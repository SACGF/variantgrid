from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import cached_property
from typing import Set, Union, Optional, Iterator

import itertools
from django.contrib.auth.models import User
from django.db.models import Q, QuerySet
from django.http import HttpRequest
from more_itertools.more import peekable
from classification.enums import ShareLevel, AlleleOriginBucket
from classification.models import ClassificationGrouping, ImportedAlleleInfo, \
    OverlapContribution
from classification.models import EvidenceKeyMap
from genes.hgvs import CHGVS
from library.utils import local_date_string
from snpdb.models import Organization, Lab, GenomeBuild, Variant, Allele
import re


@dataclass
class ClassificationGroupingExportFilter:
    user: User
    lab_mode: str
    exclude_sources: Optional[Set[Union[Lab, Organization]]] = None
    include_sources: Optional[Set[Lab]] = None
    since: Optional[datetime] = None  # has to work on classification grouping and conflict lab
    allele_origin: Optional[AlleleOriginBucket] = None

    @cached_property
    def date_str(self):
        return local_date_string()

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationGroupingExportFilter':
        exclude_sources = None
        include_sources = None
        allele_origin = None
        lab_mode = "all"
        if labs_str := request.GET.get("labs"):
            lab_mode = labs_str
            match labs_str:
                case "exclude-my-labs":
                    exclude_sources = set(Lab.valid_labs_qs(request.user, admin_check=False))
                case "include-only-my-labs":
                    include_sources = set(Lab.valid_labs_qs(request.user, admin_check=False))
                case "include-all-labs":
                    pass
                case _:
                    lab_mode = "custom-labs"
                    diff_labs = labs_str.split(",")
                    include_sources = Lab.objects.filter(group_name__in=diff_labs)

        # genome_build = GenomeBuild.get_name_or_alias(request.GET.get("genome_build"))

        since = None
        if since_str := request.GET.get("since"):
            if re.match("^[0-9]+$", since_str):
                since_days = int(since_str)
                # TODO round off to midnight
                since = datetime.now() - timedelta(days=since_days)
                since = since.replace(hour=0, minute=0, second=0, microsecond=0)
            elif date_match := re.match(r"(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})", since_str):
                since = datetime(year=int(date_match.group("year")), month=int(date_match.group("month")), day=int(date_match.group("day")))

        if allele_origin_str := request.GET.get("allele_origin"):
            allele_origin = AlleleOriginBucket(allele_origin_str)

        return ClassificationGroupingExportFilter(
            user=request.user,
            lab_mode=lab_mode,
            exclude_sources=exclude_sources,
            include_sources=include_sources,
            since=since,
            allele_origin=allele_origin
        )

    def queryset(self, genome_build: Optional[GenomeBuild] = None) -> QuerySet[ClassificationGrouping]:
        groupings = ClassificationGrouping.objects.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, classification_count__gt=0)
        if exclude_sources := self.exclude_sources:
            if exclude_orgs := [item for item in exclude_sources if isinstance(item, Organization)]:
                groupings = groupings.filter(lab__organization__in=exclude_orgs)
            if exclude_labs := [item for item in exclude_sources if isinstance(item, Lab)]:
                groupings = groupings.filter(lab__in=exclude_labs)
        if include_sources := self.include_sources:
            groupings = groupings.filter(lab__in=include_sources)

        if since := self.since:
            # update based on either ConflictLab updating, or the latest classification of a grouping updating
            # and the grouping itself updating (make sure we don't do any null updates)

            via_updated_conflicts = OverlapContribution.objects.filter(modified__gte=since).values_list("classification_grouping", flat=True)
            groupings = groupings.filter(Q(latest_classification_modification__modified__gte=since) | Q(pk__in=via_updated_conflicts))

        # order by allele ordering
        if not genome_build:
            genome_build = GenomeBuild.grch38()

        if allele_origin := self.allele_origin:
            groupings = groupings.filter(allele_origin_grouping__allele_origin_bucket=allele_origin)

        allele_sort_column = ImportedAlleleInfo.column_name_for_build(genome_build, "latest_allele_info", "genomic_sort")
        # could we build the ordering into AlleleOriginGrouping?
        groupings = groupings.order_by(
            allele_sort_column,
            'allele_origin_grouping__pk',
            'allele_origin_grouping__allele_origin_bucket',
            'allele_origin_grouping__testing_context_bucket',
            'allele_origin_grouping__tumor_type_category')

        groupings = groupings.select_related(
            "allele_origin_grouping",
            "allele_origin_grouping__allele",
            "latest_classification_modification__classification"
        )

        return groupings

    def latest_date(self) -> Optional[datetime]:
        if latest_record := self.queryset().order_by('-latest_classification_modification__modified').values_list('latest_classification_modification__modified', flat=True).first():
            return latest_record
        return None


@dataclass
class ClassificationGroupingExportFileSettings:
    rows_per_file: Optional[int] = None
    file_prefix: str = "classification_groups"
    file_include_date: bool = True

    @classmethod
    def from_request(cls, request) -> 'ClassificationGroupingExportFileSettings':
        rows_per_file: Optional[int] = None
        if rows_per_file_str := request.GET.get("rows_per_file"):
            rows_per_file = int(rows_per_file_str)
        return ClassificationGroupingExportFileSettings(
            rows_per_file=rows_per_file
        )


@dataclass
class ClassificationGroupingExportFormatProperties:
    http_content_type: str = "text/html"
    extension: str = "txt"
    delimiter_for_row: str = "\n"


@dataclass(frozen=True)
class ClassificationGroupingByAllele:
    allele_id: int
    classification_groupings: list[ClassificationGrouping]
    genome_build: Optional[GenomeBuild]

    @cached_property
    def variant(self) -> Optional[Variant]:
        if genome_build := self.genome_build:
            allele = Allele.objects.get(pk=self.allele_id)
            variant = allele.variant_for_build_optional(genome_build)
            return variant
        return None

    @cached_property
    def c_hgvs(self) -> Optional[CHGVS]:
        for cg in self.classification_groupings:
            if allele_info := cg.latest_allele_info:
                if preferred_build := allele_info[self.genome_build]:
                    if c_hgvs_obj := preferred_build.c_hgvs_obj:
                        return c_hgvs_obj
        return None

    def sub_by_allele_origin(self) -> list['ClassificationGroupingByAlleleAndOrigin']:
        by_allele_origin = defaultdict(list)
        for grouping in self.classification_groupings:
            by_allele_origin[grouping.allele_origin_bucket].append(grouping)
        return [ClassificationGroupingByAlleleAndOrigin(self, allele_origin, groupings) for allele_origin, groupings in by_allele_origin.items()]


@dataclass(frozen=True)
class ClassificationGroupingByAlleleAndOrigin:
    grouping_by_allele: ClassificationGroupingByAllele
    allele_origin_bucket: AlleleOriginBucket
    classification_groupings: list[ClassificationGrouping]

    @property
    def allele_id(self):
        return self.grouping_by_allele.allele_id

    @property
    def variant(self) -> Optional[Variant]:
        return self.grouping_by_allele.variant

    @property
    def genome_build(self) -> GenomeBuild:
        return self.grouping_by_allele.genome_build

    @property
    def c_hgvs(self):
        return self.grouping_by_allele.c_hgvs


class ClassificationGroupingExportFormat(ABC):

    def __init__(self, classification_grouping_filter: ClassificationGroupingExportFilter):
        self.classification_grouping_filter = classification_grouping_filter

    @classmethod
    @abstractmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        raise NotImplementedError()

    @property
    def genome_build(self) -> 'GenomeBuild':
        # implement if format properties state that genome build is important
        raise NotImplementedError()

    @cached_property
    def e_keys(self) -> 'EvidenceKeyMap':
        return EvidenceKeyMap.cached()

    def queryset(self, genome_build: Optional[GenomeBuild] = None) -> QuerySet[ClassificationGrouping]:
        return self.classification_grouping_filter.queryset(genome_build=genome_build)

    def allele_group_iterator(self) -> Iterator[ClassificationGroupingByAllele]:
        for allele_id, cgs in itertools.groupby(self.queryset(self.genome_build).iterator(), lambda cg: cg.allele_origin_grouping.allele.pk):
            yield ClassificationGroupingByAllele(
                allele_id,
                list(cgs),
                genome_build=self.genome_build
            )

    def header(self) -> list[str]:
        return []

    def footer(self) -> list[str]:
        return []

    @abstractmethod
    def single_row_generator(self) -> Iterator[str]:
        raise NotImplementedError()

    def row_generator(self) -> Iterator[list[str]]:
        # while we migrate from being able to return multiple rows
        for row in self.single_row_generator():
            yield [row]

    def peekable(self) -> peekable:
        return peekable(self.row_generator())

    def extra_filename_parts(self) -> list[str]:
        return []
