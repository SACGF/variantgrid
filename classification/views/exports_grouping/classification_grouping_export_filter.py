from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Set, Union, Optional, Iterator
from django.contrib.auth.models import User
from django.db.models import Q
from django.http import HttpRequest
from more_itertools.more import peekable
from classification.enums import ShareLevel
from classification.models import ClassificationGrouping, ConflictLab, ImportedAlleleInfo
from snpdb.models import Organization, Lab, GenomeBuild
import re


@dataclass
class ClassificationGroupingExportFilter:
    user: User
    genome_build: GenomeBuild
    exclude_sources: Optional[Set[Union[Lab, Organization]]] = None
    include_sources: Optional[Set[Lab]] = None
    since: Optional[datetime] = None  # has to work on classification grouping and conflict lab

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationGroupingExportFilter':
        exclude_sources = None
        include_sources = None
        if labs_str := request.GET.get("labs"):
            match labs_str:
                case "exclude-my-labs":
                    exclude_sources = set(Lab.valid_labs_qs(request.user, admin_check=False))
                case "include-only-my-labs":
                    include_sources = set(Lab.valid_labs_qs(request.user, admin_check=False))
                case "include-all-labs":
                    pass
                case _:
                    diff_labs = labs_str.split(",")
                    include_sources = Lab.objects.filter(group_name__in=diff_labs)

        genome_build = GenomeBuild.get_name_or_alias(request.GET.get("genome_build"))

        since = None
        if since_str := request.GET.get("since"):
            if re.match("[0-9]+", since_str):
                since_days = int(since_str)
                # TODO round off to midnight
                since = datetime.now() - timedelta(days=since_days)
            elif date_match := re.match(r"(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})", since_str):
                since = datetime(year=int(date_match.group("year")), month=int(date_match.group("months")), day=int(date_match.group("days")))
                # TODO truncate date

        return ClassificationGroupingExportFilter(
            user=request.user,
            genome_build=genome_build,
            exclude_sources=exclude_sources,
            include_sources=include_sources,
            since=since
        )

    def queryset(self):
        groupings = ClassificationGrouping.objects.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)
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

            via_updated_conflicts = ConflictLab.objects.filter(modified__gte=since).values_list("classification_grouping", flat=True)
            groupings = groupings.filter(Q(latest_classification_modification__modified__gte=since) | Q(pk__in=via_updated_conflicts))

        # order by allele ordering
        allele_sort_column = ImportedAlleleInfo.column_name_for_build(self.genome_build, "latest_allele_info", "genomic_sort")
        # could we build the ordering into AlleleOriginGrouping?
        groupings = groupings.order_by(
            allele_sort_column,
            'allele_origin_grouping__allele_origin_bucket',
            'allele_origin_grouping__testing_context_bucket',
            'allele_origin_grouping__tumor_type_category')

        groupings = groupings.select_related(
            "allele_origin_grouping",
            "latest_classification_modification__classification"
        )

        return groupings


@dataclass
class ClassificationGroupingExportFileSettings:
    rows_per_file: Optional[int] = None
    file_prefix: str = "classification_groups"
    file_include_date: bool = True
    row_limit: Optional[int] = None


@dataclass
class ClassificationGroupingExportFormatProperties:
    is_genome_build_relevant: bool = True
    filename_suffix: str = ""
    http_content_type: str = "text/html"
    extension: str = "txt"
    delimiter_for_row: str = "\n"


class ClassificationGroupingExportFormat(ABC):

    def __init__(self, classification_grouping_filter: ClassificationGroupingExportFilter):
        self.classification_grouping_filter = classification_grouping_filter

    @classmethod
    @abstractmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        raise NotImplementedError()

    def queryset(self):
        return self.classification_grouping_filter.queryset()

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

