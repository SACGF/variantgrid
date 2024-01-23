from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from functools import cached_property, reduce
from operator import __or__
from typing import Type, Union, Optional, Iterator, Any, Callable

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.http import HttpRequest
from guardian.shortcuts import get_objects_for_user
from threadlocals.threadlocals import get_current_request

from classification.enums import ShareLevel, ClinicalContextStatus
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import ClassificationModification, Classification, classification_flag_types, \
    DiscordanceReport, ClinicalContext, ImportedAlleleInfo
from flags.models import FlagsMixin, Flag, FlagComment
from library.utils import batch_iterator, local_date_string, http_header_date_now
from snpdb.models import GenomeBuild, Lab, Organization, Allele, Variant


@dataclass
class ClassificationIssue:
    classification: ClassificationModification
    withdrawn = False
    validation_include = False
    not_matched = False

    @property
    def message(self) -> Optional[str]:
        messages: list[str] = []
        if self.withdrawn:
            messages.append("Classification has been withdrawn")
        # if self.transcript_version or self.matching_warning:
        #     messages.append("Requires confirmation of variant match")
        # if self.c_37_not_38:
        #     messages.append("transcript across genome builds requires confirmation")
        if not self.validation_include:
            messages.append("Excluded due to outstanding variant matching/liftover warning")
        if self.not_matched:
            messages.append("Could not liftover/normalise")

        if messages:
            return ", ".join(messages)
        else:
            return None

    def __repr__(self):
        return self.message

    @property
    def has_issue(self):
        return self.withdrawn or self.not_matched or not self.validation_include


@dataclass
class AlleleData:
    """
    All the classifications for a single Allele within an export
    :var source: The source of all export data
    :var allele_id: The Allele ID that all the classifications belong to
    :var all_cms: All the classifications that passed the filter, but might include rows with errors that most
    export will want to ignore (see cms instead)
    """
    source: 'ClassificationFilter'
    allele_id: int
    all_cms: list[ClassificationIssue] = field(default_factory=list)  # misleading name, should be all_ci or something

    def sort(self):
        def record_order(ci: ClassificationIssue):
            return ci.classification.lab.group_name, ci.classification.pk

        self.all_cms.sort(key=record_order)

    cached_allele: Optional[Allele] = None
    cached_variant: Optional[Variant] = None
    cached_data: dict[str, Any] = None

    @staticmethod
    def from_allele_info(source: 'ClassificationFilter', allele_info: ImportedAlleleInfo):
        variant: Optional[Variant] = None
        if variant_info := allele_info[source.genome_build]:
            variant = variant_info.variant
        return AlleleData(
            source=source,
            allele_id=allele_info.allele_id,
            cached_allele=allele_info.allele,
            cached_variant=variant
        )

    def __setitem__(self, key, value):
        if not self.cached_data:
            self.cached_data = {}
        self.cached_data[key] = value

    def __getitem__(self, item):
        if cached_data := self.cached_data:
            return cached_data.get(item)

    def __contains__(self, item):
        if cached_data := self.cached_data:
            return item in cached_data
        return False

    # @staticmethod
    # def pre_process(batch: list['AlleleData']):
    #     if batch:
    #         genome_build = batch[0].genome_build
    #         allele_ids = [ad.allele_id for ad in batch]
    #         alleles = Allele.objects.filter(pk__in=allele_ids).select_related('clingen_allele')
    #         allele_dict = {allele.pk: allele for allele in alleles}
    #         # FIXME select related all this variant stuff outside?
    #         variant_alleles = VariantAllele.objects.filter(allele_id__in=allele_ids, genome_build=genome_build).select_related('variant', 'variant__locus', 'variant__locus__contig', 'variant__locus__ref', 'variant__alt')
    #         variant_dict = {variant_allele.allele_id: variant_allele.variant for variant_allele in variant_alleles}
    #
    #         for ad in batch:
    #             ad.cached_allele = allele_dict.get(ad.allele_id)
    #             ad.cached_variant = variant_dict.get(ad.allele_id)

    @property
    def allele(self) -> Allele:
        # return Allele.objects.filter(pk=self.allele_id).select_related('clingen_allele').first()
        return self.cached_allele

    @cached_property
    def variant(self) -> Optional[Variant]:
        return self.cached_variant
        # if allele := self.allele:
        #    try:
        #        return self.allele.variant_for_build(genome_build=self.genome_build, best_attempt=True)
        #    except ValueError:
        #        pass

    @property
    def genome_build(self) -> GenomeBuild:
        return self.source.genome_build

    @property
    def cms(self) -> list[ClassificationModification]:
        # The classifications that should be exported (passed validation, not withdrawn)
        return [ci.classification for ci in self.all_cms if not ci.has_issue]

    @property
    def cms_regardless_of_issues(self) -> list[ClassificationModification]:
        return [ci.classification for ci in self.all_cms]

    @property
    def issues(self) -> list[ClassificationIssue]:
        # All Classification Issues that actually have issues
        return [ci for ci in self.all_cms if ci.has_issue]

    def __bool__(self):
        return bool(self.all_cms)


def flag_ids_to(model: Type[FlagsMixin], qs: Union[QuerySet[Flag], QuerySet[FlagComment]]) -> set[int]:
    """
    Convert a qs of Flags or FlagComments to ids of Classification/Allele/etc.
    :param model: The Model to get the IDs of e.g. Classification, Allele, assumes the flag qs was for Flags that
    correspond to the model.
    :param qs: A QuerySet of relevant flags e.g. flags_ids_to(model: Classification, qs: ...open transcript mismatch flag...)
    :return: A set of model IDs that are linked to the flag qs
    """
    if qs.model == Flag:
        qs = qs.values_list('collection_id', flat=True)\
            .order_by('collection_id').distinct('collection_id')
    elif qs.model == FlagComment:
        qs = qs.values_list('flag__collection_id', flat=True)\
            .order_by('flag__collection_id').distinct('flag__collection_id')
    else:
        raise ValueError(f"Cannot deal with a QuerySet of model {qs.model}")

    return set(model.objects.filter(flag_collection_id__in=qs).values_list('pk', flat=True))


class TranscriptStrategy(str, Enum):
    """
    What transcripts are going to be included
    :cvar ALL: Include all transcripts
    :cvar REFSEQ: More accurately only include transcripts that Alissa supports NM/NR
    """
    ALL = "all"
    REFSEQ = "refseq"


class DiscordanceReportStatus(str, Enum):
    """
    Different from DiscordanceReportResolution as that uses None to indicate an ongoing discordance
    where when we need it that would indicate no discordance at all
    """
    ON_GOING = 'ongoing'
    CONTINUED = 'continued'
    PENDING_CONCORDANCE = 'pending'


@dataclass
class ClassificationFilter:
    """
    Attributes:
        user: The user performing the export
        genome_build: The genome build the user has requested (not always relevant)
        exclude_sources: Labs/Orgs that should be excluded
        include_sources: If provided Labs/Orgs that should only be included
        since: Try to optimise data so we only include alleles that have been modified since this date
        min_share_level: Minimum share level of classifications to include
        transcript_strategy: Which classification transcripts should be considered
        rows_per_file: Max number of rows per file (data from the 1 allele will still be grouped together if formatted produces multiple lines)
        starting_query: If this is provided then instead of using all of the above filters, just this
        TODO rename row_limit to indicate that it's the max per file, not
    """

    user: User
    genome_build: GenomeBuild
    exclude_sources: Optional[set[Union[Lab, Organization]]] = None
    include_sources: Optional[set[Lab]] = None
    since: Optional[datetime] = None
    min_share_level: ShareLevel = ShareLevel.LAB
    transcript_strategy: TranscriptStrategy = TranscriptStrategy.ALL
    rows_per_file: Optional[int] = None
    allele: Optional[int] = None
    file_prefix: str = "classifications"
    file_include_date: bool = True
    starting_query: Optional[QuerySet[ClassificationModification]] = None
    benchmarking: bool = False
    path_info: Optional[str] = None
    request_params: Optional[dict] = None
    row_limit: Optional[int] = None
    _last_modified: str = None

    def __post_init__(self):
        self._last_modified = http_header_date_now()
        if self.path_info and self.request_params:
            pass
        elif request := get_current_request():
            self.record_request_details(request)

    @property
    def last_modified_header(self) -> str:
        return self._last_modified

    def record_request_details(self, request):
        if not self.path_info:
            self.path_info = request.path_info
        if self.request_params is None:
            request_params = {}
            if GET := request.GET:
                for key, value in GET.items():
                    request_params[key] = value
            if POST := request.POST:
                for key, value in POST.items():
                    request_params[key] = value
            self.request_params = request_params

    @cached_property
    def date_str(self) -> str:
        return local_date_string()

    @staticmethod
    def _string_to_group_name(model: Type[Union[Lab, Organization]], group_names: str) -> Union[set[Lab], set[Organization]]:
        """
        Converts comma separated group names to the org or lab objects.
        Invalid org/lab group names will be ignored
        :param model: Organization or Lab
        :param group_names: e.g. "org_1, org_2" or "org_1/lab_x, org_1/lab_y"
        :return: A set of models that matches the group names
        """
        if not group_names:
            return set()
        parts = [gn.strip() for gn in group_names.split(',')]
        objs = [model.objects.filter(group_name=group_name).first() for group_name in parts]
        objs = [m for m in objs if m]
        return set(objs)

    @staticmethod
    def from_request(request: HttpRequest) -> 'ClassificationFilter':
        """
        Create a filter of Classification data from classification_export.html
        :param request: Request
        :return: Populated ClassificationFilter
        """
        user = request.user

        exclude_sources = \
            ClassificationFilter._string_to_group_name(Organization, request.query_params.get('exclude_orgs')) | \
            ClassificationFilter._string_to_group_name(Lab, request.query_params.get('exclude_labs'))
        include_sources = ClassificationFilter._string_to_group_name(Lab, request.query_params.get('include_labs'))
        build_name = request.query_params.get('build', 'GRCh38')

        share_level_str = request.query_params.get('share_level', 'logged_in_users')
        # public is ever so slowly being deprecated
        if share_level_str == 'public':
            share_level_str = 'logged_in_users'
        elif share_level_str == 'any':
            share_level_str = 'lab'
        share_level = ShareLevel(share_level_str)
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        transcript_strategy = TranscriptStrategy(request.query_params.get('transcript_strategy', 'all'))
        since: Optional[datetime] = None
        if since_str := request.query_params.get('since'):
            from classification.views.classification_export_view import parse_since
            since = parse_since(since_str)

        # just for debugging purposes allows us to download a single allele
        allele: Optional[int] = None
        if allele_str := request.query_params.get('allele'):
            allele = int(allele_str)

        benchmarking = request.query_params.get('benchmark') == 'true'

        rows_per_file = None
        if rows_per_file_str := request.query_params.get('rows_per_file'):
            try:
                rows_per_file = int(rows_per_file_str)
                if rows_per_file <= 0:
                    raise ValueError("Can't handle negative rows per file")
                if 0 < rows_per_file < 100:
                    rows_per_file = 100
            except:
                pass
        elif request.query_params.get("type") == "mvl":
            rows_per_file = 10000

        row_limit = None
        if row_limit_str := request.query_params.get('row_limit'):
            row_limit = int(row_limit_str)

        return ClassificationFilter(
            user=user,
            exclude_sources=exclude_sources,
            include_sources=include_sources,
            genome_build=genome_build,
            min_share_level=share_level,
            transcript_strategy=transcript_strategy,
            since=since,
            allele=allele,
            benchmarking=benchmarking,
            rows_per_file=rows_per_file,
            row_limit=row_limit
        )

    @cached_property
    def c_hgvs_col(self):
        """
        The classification column that represents the genome build that the user requested
        :return: The name of the c.hgvs column to use in a Classification QS
        """
        if self.genome_build == GenomeBuild.grch37():
            return 'classification__allele_info__grch37__c_hgvs'
        else:
            return 'classification__allele_info__grch38__c_hgvs'

    # @cached_property
    # def _bad_allele_transcripts(self) -> dict[int, set[str]]:
    #     """
    #     :return: A dictionary of Allele ID to a set of Transcripts for that allele which has bad flags
    #     """
    #
    #     qs = Flag.objects.filter(
    #         flag_type=allele_flag_types.allele_37_not_38,
    #         resolution__status=FlagStatus.OPEN
    #     ).values_list('collection_id', 'data')
    #
    #     allele_to_bad_transcripts: dict[int, set[str]] = defaultdict(set)
    #     collection_id_to_transcript = defaultdict(set)
    #     for collection_id, data in qs:
    #         if data:
    #             collection_id_to_transcript[collection_id].add(data.get('transcript'))
    #
    #     allele_qs = Allele.objects.filter(flag_collection__in=collection_id_to_transcript.keys())\
    #         .values_list('pk', 'flag_collection_id')
    #     for pk, collection_id in allele_qs:
    #         if transcripts := collection_id_to_transcript.get(collection_id):
    #             allele_to_bad_transcripts[pk].update(transcripts)
    #
    #     return allele_to_bad_transcripts
    #
    # @cached_property
    # def _transcript_version_classification_ids(self) -> set[int]:
    #     """
    #     Returns a set of classification IDs that have flags that make us want to exclude due to errors
    #     :return: A set of classification IDs
    #     """
    #     return flag_ids_to(
    #         Classification,
    #         Flag.objects.filter(
    #             flag_type=classification_flag_types.transcript_version_change_flag,
    #             resolution__status=FlagStatus.OPEN
    #         ))
    #
    # @cached_property
    # def _variant_matching_classification_ids(self) -> set[int]:
    #     return flag_ids_to(
    #         Classification,
    #         Flag.objects.filter(
    #             flag_type=classification_flag_types.matching_variant_warning_flag,
    #             resolution__status=FlagStatus.OPEN
    #         ))

    @cached_property
    def _discordant_classification_ids(self) -> dict[int, DiscordanceReportStatus]:
        """
        Returns a set of classification IDs where the classification id discordant
        Ids are not necessarily part of this import
        :return: A set of classification IDs
        """
        discordance_status: dict[int, DiscordanceReportStatus] = {}
        for cc in ClinicalContext.objects.filter(status=ClinicalContextStatus.DISCORDANT):
            if dr := DiscordanceReport.latest_report(cc):
                status: Optional[DiscordanceReportStatus]
                if dr.is_pending_concordance:
                    status = DiscordanceReportStatus.PENDING_CONCORDANCE
                elif dr.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
                    status = DiscordanceReportStatus.CONTINUED
                else:
                    status = DiscordanceReportStatus.ON_GOING

                for c in dr.actively_discordant_classification_ids():
                    discordance_status[c] = status

        return discordance_status

    def is_discordant(self, cm: ClassificationModification) -> DiscordanceReportStatus:
        if settings.DISCORDANCE_ENABLED:
            return self._discordant_classification_ids.get(cm.classification_id)

    @cached_property
    def _since_flagged_classification_ids(self) -> set[int]:
        """
        TODO rename to indicate this is a flag check only
        Returns classification ids that have had flags change since the since date
        :return: A set of classification Ids that have relevant flags that have changed since the since date
        """
        return flag_ids_to(Classification, FlagComment.objects.filter(
            flag__flag_type__in={
                classification_flag_types.classification_withdrawn,
                # classification_flag_types.transcript_version_change_flag,
                # classification_flag_types.matching_variant_warning_flag,
                classification_flag_types.classification_pending_changes
            },
            created__gte=self.since
        ))

    def _passes_since(self, allele_data: AlleleData) -> bool:
        """
        Is there anything about this AlleleData that indicates it should be included since the since date
        """
        if not self.since:
            return True
        for cmi in allele_data.all_cms:
            cm = cmi.classification
            if cm.modified >= self.since or cm.classification.modified > self.since:
                return True
            if cm.classification_id in self._since_flagged_classification_ids:
                return True
            if cm.classification.allele_info.latest_validation.modified > self.since:
                return True
        return False

    @property
    def _share_levels(self) -> set[ShareLevel]:
        """
        :return: A set of all ShareLevels that should be considered
        """
        share_levels: set[ShareLevel] = set()
        for sl in ShareLevel.ALL_LEVELS:
            if sl >= self.min_share_level:
                share_levels.add(sl)
        return share_levels

    @cached_property
    def cms_qs(self) -> QuerySet[ClassificationModification]:
        """
        Returns a new QuerySet of all classifications BEFORE
        filtering for errors and since date
        """
        if starting_query := self.starting_query:
            cms = starting_query
        else:
            cms = ClassificationModification.objects.filter(is_last_published=True)

        if self.min_share_level != ShareLevel.LAB:
            cms = cms.filter(share_level__in=self._share_levels)

        if not self.since:
            # only worry about withdrawn if doing 'since' (as we might need to report the withdrawing (json),
            # or at least be aware of it for changes (mvl))
            cms = cms.exclude(classification__withdrawn=True)

        if labs := self.include_sources:
            cms = cms.filter(classification__lab__in=labs)
        elif self.exclude_sources:
            if exclude_orgs := [org for org in self.exclude_sources if isinstance(org, Organization)]:
                cms = cms.exclude(classification__lab__organization__in=exclude_orgs)
            if exclude_labs := [lab for lab in self.exclude_sources if isinstance(lab, Lab)]:
                cms = cms.exclude(classification__lab__in=exclude_labs)

        genome_build_str = '37'
        if self.genome_build.is_version(38):
            genome_build_str = '38'

        genomic_sort = f'classification__allele_info__grch{genome_build_str}__genomic_sort'
        cms = cms.order_by(genomic_sort, 'classification__allele_info__allele', 'classification__lab', '-classification__id')

        cms = cms.select_related(
            'classification',
            'classification__allele_info__latest_validation',
            'classification__allele_info__allele',
            'classification__allele_info__allele__clingen_allele',
            f'classification__allele_info__grch{genome_build_str}__variant',
            f'classification__allele_info__grch{genome_build_str}__genome_build',
            f'classification__allele_info__grch{genome_build_str}__variant__locus',
            f'classification__allele_info__grch{genome_build_str}__variant__locus__contig',
            f'classification__allele_info__grch{genome_build_str}__variant__locus__ref',
            f'classification__allele_info__grch{genome_build_str}__variant__alt',
            'classification__lab__organization',
            'classification__clinical_context'
        )

        if allele_id := self.allele:
            cms = cms.filter(classification__allele_info__allele_id=allele_id)

        # PERMISSION CHECK
        cms = get_objects_for_user(self.user, ClassificationModification.get_read_perm(), cms, accept_global_perms=True)
        # can't filter out transcript versions (due to 37!=38) easily here, do that later

        if self.transcript_strategy == TranscriptStrategy.REFSEQ:
            from classification.views.classification_export_view import ALISSA_ACCEPTED_TRANSCRIPTS
            acceptable_transcripts: list[Q] = [
                Q(**{f'{self.c_hgvs_col}__startswith': tran}) for tran in ALISSA_ACCEPTED_TRANSCRIPTS
            ]
            cms = cms.filter(reduce(__or__, acceptable_transcripts))

        return cms

    def _record_issues(self, allele_id: int, cm: ClassificationModification) -> ClassificationIssue:
        ci = ClassificationIssue(classification=cm)
        ci.withdrawn = cm.classification.withdrawn
        # ci.transcript_version = cm.classification_id in self._transcript_version_classification_ids
        # ci.matching_warning = cm.classification_id in self._variant_matching_classification_ids
        if not allele_id or not ci.classification.classification.get_c_hgvs(self.genome_build):
            ci.not_matched = True
        if (allele_info := cm.classification.allele_info) and (latest_validation := allele_info.latest_validation):
            ci.validation_include = latest_validation.include

        return ci

    def _allele_data(self) -> Iterator[AlleleData]:
        """
        Convert the classification query set into AlleleData
        Is not filtered at this point
        """
        allele_data: Optional[AlleleData] = None
        for cm in self.cms_qs.iterator(chunk_size=1000):
            if allele_info := cm.classification.allele_info:
                allele_id = cm.classification.allele_id

                if not allele_data or allele_id != allele_data.allele_id:
                    if allele_data:
                        allele_data.sort()
                        yield allele_data

                    allele_data = AlleleData.from_allele_info(source=self, allele_info=cm.classification.allele_info)
                allele_data.all_cms.append(self._record_issues(allele_id=allele_id, cm=cm))

        if allele_data:
            allele_data.sort()
            yield allele_data

    def _allele_data_filtered(self) -> Iterator[AlleleData]:
        """
        Main method of getting data out of ClassificationFilter
        Is only AlleleData with classifications that have passed since checks and errors
        """
        for allele_data in self._allele_data():
            if self._passes_since(allele_data):
                yield allele_data

    def allele_data_filtered_pre_processed(self, batch_processor: Optional[Callable[[list[AlleleData]], None]] = None) -> Iterator[AlleleData]:
        for batch in batch_iterator(self._allele_data_filtered(), batch_size=100):
            # AlleleData.pre_process(batch)
            if batch_processor:
                batch_processor(batch)
            for data in batch:
                yield data
