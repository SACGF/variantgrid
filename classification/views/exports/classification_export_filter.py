from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from functools import reduce
from operator import __or__
from typing import List, Type, Union, Set, Optional, Dict, Iterator
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.http import HttpRequest
from django.utils.timezone import now
from guardian.shortcuts import get_objects_for_user
from lazy import lazy
from pytz import timezone

from classification.enums import ShareLevel
from classification.models import ClassificationModification, Classification, classification_flag_types
from flags.models import FlagsMixin, Flag, FlagComment, FlagStatus
from snpdb.models import GenomeBuild, Lab, Organization, allele_flag_types, Allele


@dataclass
class ClassificationIssue:
    classification: ClassificationModification
    withdrawn = False
    transcript_version = False
    matching_warning = False
    c_37_not_38 = False
    not_matched = False

    @property
    def message(self) -> Optional[str]:
        messages: List[str] = list()
        if self.withdrawn:
            messages.append("Classification has been withdrawn")
        if self.transcript_version or self.matching_warning:
            messages.append("Requires confirmation of variant match")
        if self.c_37_not_38:
            messages.append("transcript across genome builds requires confirmation")
        if self.not_matched:
            messages.append("Could not liftover/normalise")
        if messages:
            return ", ".join(messages)
        else:
            return None

    @property
    def has_issue(self):
        return self.withdrawn or self.transcript_version or self.matching_warning or self.c_37_not_38 or self.not_matched


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
    all_cms: List[ClassificationIssue] = field(default_factory=list)

    @property
    def cms(self) -> List[ClassificationModification]:
        # The classifications that should be exported (passed validation, not withdrawn)
        return [ci.classification for ci in self.all_cms if not ci.has_issue]

    @property
    def issues(self) -> List[ClassificationIssue]:
        # All Classification Issues that actually have issues
        return [ci for ci in self.all_cms if ci.has_issue]

    def __bool__(self):
        return bool(self.all_cms)


def flag_ids_to(model: Type[FlagsMixin], qs: Union[QuerySet[Flag], QuerySet[FlagComment]]) -> Set[int]:
    """
    Convert a qs of Flags or FlagComments to ids of Classification/Allele/etc
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
        row_limit: Max number of rows per file (data from the 1 allele will still be grouped together if formatted produces multiple lines)
        starting_query: If this is provided then instead of using all of the above filters, just this
        TODO rename row_limit to indicate that it's the max per file, not
    """

    user: User
    genome_build: GenomeBuild
    exclude_sources: Optional[Set[Union[Lab, Organization]]] = None
    include_sources: Optional[Set[Lab]] = None
    since: Optional[datetime] = None
    min_share_level: ShareLevel = ShareLevel.LAB
    transcript_strategy: TranscriptStrategy = TranscriptStrategy.ALL
    rows_per_file: Optional[int] = None
    allele: Optional[int] = None
    file_prefix: str = "classifications"
    file_include_date: bool = True
    starting_query: Optional[QuerySet[Classification]] = None
    benchmarking: bool = False

    @lazy
    def date_str(self) -> str:
        return now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")

    @staticmethod
    def _string_to_group_name(model: Type[Union[Lab, Organization]], group_names: str) -> Union[Set[Lab], Set[Organization]]:
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

        # TODO include rows_per_file into filter? right now it's hardcoded when doing MVL

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
            rows_per_file=rows_per_file
        )

    @lazy
    def c_hgvs_col(self):
        """
        The classification column that represents the genome build that the user requested
        :return: The name of the c.hgvs column to use in a Classification QS
        """
        if self.genome_build == GenomeBuild.grch37():
            return 'classification__chgvs_grch37'
        else:
            return 'classification__chgvs_grch38'

    @lazy
    def _bad_allele_transcripts(self) -> Dict[int, Set[str]]:
        """
        :return: A dictionary of Allele ID to a set of Transcripts for that allele which has bad flags
        """

        qs = Flag.objects.filter(
            flag_type=allele_flag_types.allele_37_not_38,
            resolution__status=FlagStatus.OPEN
        ).values_list('collection_id', 'data')

        allele_to_bad_transcripts: Dict[int, Set[str]] = defaultdict(set)
        collection_id_to_transcript = defaultdict(set)
        for collection_id, data in qs:
            if data:
                collection_id_to_transcript[collection_id].add(data.get('transcript'))

        allele_qs = Allele.objects.filter(flag_collection__in=collection_id_to_transcript.keys())\
            .values_list('pk', 'flag_collection_id')
        for pk, collection_id in allele_qs:
            if transcripts := collection_id_to_transcript.get(collection_id):
                allele_to_bad_transcripts[pk].update(transcripts)

        return allele_to_bad_transcripts

    @lazy
    def _transcript_version_classification_ids(self) -> Set[int]:
        """
        Returns a set of classification IDs that have flags that make us want to exclude due to errors
        :return: A set of classification IDs
        """
        return flag_ids_to(
            Classification,
            Flag.objects.filter(
                flag_type=classification_flag_types.transcript_version_change_flag,
                resolution__status=FlagStatus.OPEN
            ))

    @lazy
    def _variant_matching_classification_ids(self) -> Set[int]:
        return flag_ids_to(
            Classification,
            Flag.objects.filter(
                flag_type=classification_flag_types.matching_variant_warning_flag,
                resolution__status=FlagStatus.OPEN
            ))

    @lazy
    def _discordant_classification_ids(self) -> Set[int]:
        """
        Returns a set of classification IDs where the classification id discordant
        Ids are not necessarily part of this import
        :return: A set of classification IDs
        """
        return flag_ids_to(Classification, Flag.objects.filter(
            flag_type=classification_flag_types.discordant,
            resolution__status=FlagStatus.OPEN
        ))

    def is_discordant(self, cm: ClassificationModification) -> bool:
        return cm.classification_id in self._discordant_classification_ids

    @lazy
    def _since_flagged_classification_ids(self) -> Set[int]:
        """
        TODO rename to indicate this is a flag check only
        Returns classification ids that have had flags change since the since date
        :return: A set of classification Ids that have relevant flags that have changed since the since date
        """
        return flag_ids_to(Classification, FlagComment.objects.filter(
            flag__flag_type__in={
                classification_flag_types.classification_withdrawn,
                classification_flag_types.transcript_version_change_flag,
                classification_flag_types.matching_variant_warning_flag
            },
            created__gte=self.since
        ))

    @lazy
    def _since_flagged_allele_ids(self) -> Set[int]:
        """
        :return: A set of allele IDs that have relevant flags that have changed since the since date
        """
        return flag_ids_to(Allele, FlagComment.objects.filter(
            flag__flag_type__in={
                allele_flag_types.allele_37_not_38
            },
            created__gte=self.since
        ))

    def _passes_since(self, allele_data: AlleleData) -> bool:
        """
        Is there anything about this AlleleData that indicates it should be included since the since date
        """
        if not self.since:
            return True
        if allele_data.allele_id and allele_data.allele_id in self._since_flagged_allele_ids:
            return True
        for cmi in allele_data.all_cms:
            cm = cmi.classification
            if cm.modified >= self.since or cm.classification.modified > self.since:
                return True
            if cm.classification_id in self._since_flagged_classification_ids:
                return True
        return False

    @property
    def _share_levels(self) -> Set[ShareLevel]:
        """
        :return: A set of all ShareLevels that should be considered
        """
        share_levels: Set[ShareLevel] = set()
        for sl in ShareLevel.ALL_LEVELS:
            if sl >= self.min_share_level:
                share_levels.add(sl)
        return share_levels

    def cms_qs(self) -> QuerySet[ClassificationModification]:
        """
        Returns a new QuerySet of all classifications BEFORE
        filtering for errors and since date
        """
        if starting_query := self.starting_query:
            cms = starting_query
        else:
            cms = ClassificationModification.objects.filter(is_last_published=True)

        cms = cms.filter(share_level__in=self._share_levels)

        if not self.since:
            # can only exclude withdrawn from complete consideration
            cms = cms.exclude(classification__withdrawn=True)

        if labs := self.include_sources:
            cms = cms.filter(classification__lab__in=labs)
        elif self.exclude_sources:
            if exclude_orgs := [org for org in self.exclude_sources if isinstance(org, Organization)]:
                cms = cms.exclude(classification__lab__organization__in=exclude_orgs)
            if exclude_labs := [lab for lab in self.exclude_sources if isinstance(lab, Lab)]:
                cms = cms.exclude(classification__lab__in=exclude_labs)

        cms = cms.order_by('-classification__allele_id', '-classification__id')
        cms = cms.select_related('classification', 'classification__lab', 'classification__lab__organization', 'classification__clinical_context')

        if allele_id := self.allele:
            cms = cms.filter(classification__allele_id=allele_id)

        # Always safe to exclude these (unless we want them in a CSV) even with since changes
        # couldn't show these ones if we wanted to

        ## Let these bad records in just so they can be put into errors
        ##
        # cms = cms.exclude(classification__allele__isnull=True).exclude(classification__variant__isnull=True)
        # cms = cms.exclude(**{f'{self.c_hgvs_col}__isnull': True})

        # PERMISSION CHECK
        cms = get_objects_for_user(self.user, ClassificationModification.get_read_perm(), cms, accept_global_perms=True)
        # can't filter out transcript versions (due to 37!=38) easily here, do that later

        if self.transcript_strategy == TranscriptStrategy.REFSEQ:
            from classification.views.classification_export_view import ALISSA_ACCEPTED_TRANSCRIPTS
            acceptable_transcripts: List[Q] = [
                Q(**{f'{self.c_hgvs_col}__startswith': tran}) for tran in ALISSA_ACCEPTED_TRANSCRIPTS
            ]
            cms = cms.filter(reduce(__or__, acceptable_transcripts))

        return cms

    def _record_issues(self, allele_id: int, cm: ClassificationModification) -> ClassificationIssue:
        ci = ClassificationIssue(classification=cm)
        ci.withdrawn = cm.classification.withdrawn
        ci.transcript_version = cm.classification_id in self._transcript_version_classification_ids
        ci.matching_warning = cm.classification_id in self._variant_matching_classification_ids
        if not allele_id:
            ci.not_matched = True
        else:
            if allele_bad_transcripts := self._bad_allele_transcripts.get(allele_id):
                ci.c_37_not_38 = cm.transcript in allele_bad_transcripts

        if not ci.classification.classification.get_c_hgvs(self.genome_build):
            ci.not_matched = True

        return ci

    def _allele_data(self) -> Iterator[AlleleData]:
        """
        Convert the classification query set into AlleleData
        Is not filtered at this point
        """
        allele_data: Optional[AlleleData] = None
        for cm in self.cms_qs():
            allele_id = cm.classification.allele_id
            # FIXME: make an AlleleData for no Allele
            if not allele_data or allele_id != allele_data.allele_id:
                if allele_data:
                    yield allele_data
                allele_data = AlleleData(source=self, allele_id=allele_id)
            allele_data.all_cms.append(self._record_issues(allele_id=allele_id, cm=cm))

        if allele_data:
            yield allele_data

    def allele_data_filtered(self) -> Iterator[AlleleData]:
        """
        Main method of getting data out of ClassificationFilter
        Is only AlleleData with classifications that have passed since checks and errors
        """
        for allele_data in self._allele_data():
            if self._passes_since(allele_data):
                yield allele_data
