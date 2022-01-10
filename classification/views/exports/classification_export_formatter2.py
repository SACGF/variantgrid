import zipfile
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from io import StringIO
from typing import Optional, Set, Union, List, Iterator, Dict, Tuple, Type, Iterable
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.http.response import HttpResponseBase
from django.utils.timezone import now
from pytz import timezone
from cyvcf2.cyvcf2 import defaultdict
from django.http import HttpRequest, HttpResponse, StreamingHttpResponse
from guardian.shortcuts import get_objects_for_user
from lazy import lazy
from threadlocals.threadlocals import get_current_request

from classification.enums import ShareLevel
from classification.models import ClassificationModification, classification_flag_types, Classification
from flags.models import Flag, FlagStatus, FlagComment, FlagsMixin
from genes.hgvs import CHGVS
from library.guardian_utils import bot_group
from library.log_utils import NotificationBuilder
from snpdb.models import Lab, Organization, GenomeBuild, allele_flag_types, Allele


@dataclass
class AlleleData:
    """
    All the classifications for a single Allele within an export
    :var source: The source of all export data
    :var allele_id: The Allele ID that all the classifications belong to
    :var cms: The classifications that should be exported (passed validation, not withdrawn)
    """
    source: 'ClassificationFilter'
    allele_id: int
    cms: List[ClassificationModification] = field(default_factory=list)

    def __bool__(self):
        return bool(self.cms)


@dataclass
class CHGVSData:
    """
    A sub-division of AlleleData.
    Will create one record per unique c.hgvs string within the allele*
    (c.hgvs differing in just transcript version are still bundled together)

    :var source: The allele data record (TODO rename)
    :var chgvs: The c.hgvs with the highest found transcript version
    :var different_versions: Bool indicating if multiple transcript versions were bundled together here
    :var cms: The classifications
    """
    source: AlleleData
    chgvs: CHGVS
    different_versions: bool = False
    cms: List[ClassificationModification] = field(default_factory=list)

    @staticmethod
    def split_into_c_hgvs(
            allele_data: AlleleData,
            use_full: bool) -> List['CHGVSData']:
        """
        Break up an AlleleData into sub CHGVSDatas
        :param allele_data: The Alissa data to split
        :param use_full: Should the c.hgvs use explicit bases when optional (required by Alissa)
        :return: An array of c.hgvs based data, most often will only be 1 record
        """
        genome_build = allele_data.source.genome_build
        cms = allele_data.cms
        by_chgvs: Dict[CHGVS, List[Tuple[ClassificationModification, int]]] = defaultdict(list)
        for cm in cms:
            c_hgvs = CHGVS(cm.classification.get_c_hgvs(genome_build=genome_build, use_full=use_full))
            by_chgvs[c_hgvs.without_transcript_version].append((cm, c_hgvs.transcript_parts.version))

        sub_datas: List[CHGVSData] = list()
        for c_hgvs, versioned_cms in by_chgvs.items():
            cms = [vcm[0] for vcm in versioned_cms]
            versions = [v for v in reversed(sorted(set(vcm[1] for vcm in versioned_cms))) if v is not None]
            # always provide highest transcript version when there's multiple
            # but make sure we at least have 1 transcript version and not an LRG for example
            if versions:
                c_hgvs.transcript = c_hgvs.transcript + f'.{versions[0]}'

            sub_datas.append(CHGVSData(source=allele_data, chgvs=c_hgvs, different_versions=len(cms) > 1, cms=cms))
        return sub_datas


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
        TODO rename row_limit to indicate that it's the max per file, not
    """
    user: User
    genome_build: GenomeBuild
    exclude_sources: Optional[Set[Union[Lab, Organization]]] = None
    include_sources: Optional[Set[Lab]] = None
    since: Optional[datetime] = None
    min_share_level: ShareLevel = ShareLevel.ALL_USERS
    transcript_strategy: TranscriptStrategy = TranscriptStrategy.ALL
    row_limit: Optional[int] = None

    @lazy
    def date_str(self) -> str:
        self.user
        return now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")

    @staticmethod
    def string_to_group_name(model: Type[Union[Lab, Organization]], group_names: str) -> Union[Set[Lab], Set[Organization]]:
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
            ClassificationFilter.string_to_group_name(Organization, request.query_params.get('exclude_orgs')) | \
            ClassificationFilter.string_to_group_name(Lab, request.query_params.get('exclude_labs'))
        include_sources = ClassificationFilter.string_to_group_name(Lab, request.query_params.get('include_labs'))
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

        # FIXME include row_limit into filter? right now it's hardcoded when doing MVL

        return ClassificationFilter(
            user=user,
            exclude_sources=exclude_sources,
            include_sources=include_sources,
            genome_build=genome_build,
            min_share_level=share_level,
            transcript_strategy=transcript_strategy,
            since=since
            # row_limit=100
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
    def bad_allele_transcripts(self) -> Dict[int, Set[str]]:
        """
        :return: A dictionary of Allele ID to a set of Transcripts for that allele which has bad flags
        """

        qs = Flag.objects.filter(
            flag_type=allele_flag_types.allele_37_not_38,
            resolution__status=FlagStatus.OPEN
        ).values_list('collection_id', 'data')

        allele_to_bad_transcripts: Dict[int, Set[str]] = defaultdict(set)
        collection_id_to_transcript = defaultdict(str)
        for collection_id, data in qs:
            if data:
                collection_id_to_transcript[collection_id] = data.get('transcript')

        allele_qs = Allele.objects.filter(flag_collection__in=collection_id_to_transcript.keys())\
            .values_list('pk', 'flag_collection_id')
        for pk, collection_id in allele_qs:
            if transcript := collection_id_to_transcript.get(collection_id):
                allele_to_bad_transcripts[pk].add(transcript)

        return allele_to_bad_transcripts

    @lazy
    def bad_classification_ids(self) -> Set[int]:
        """
        Returns a set of classification IDs that have flags that make us want to exclude due to errors
        :return: A set of classification IDs
        """
        return flag_ids_to(Classification, Flag.objects.filter(
            flag_type__in={classification_flag_types.transcript_version_change_flag,
                           classification_flag_types.matching_variant_warning_flag},
            resolution__status=FlagStatus.OPEN
        ))

    @lazy
    def discordant_classification_ids(self) -> Set[int]:
        """
        Returns a set of classification IDs where the classification id discordant
        :return: A set of classification IDs
        """
        return flag_ids_to(Classification, Flag.objects.filter(
            flag_type=classification_flag_types.discordant,
            resolution__status=FlagStatus.OPEN
        ))

    def is_discordant(self, cm: ClassificationModification) -> bool:
        return cm.classification_id in self.discordant_classification_ids

    @lazy
    def since_classifications(self) -> Set[int]:
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
    def since_alleles(self) -> Set[int]:
        """
        TODO rename to indicate this is a flag check only
        :return: A set of allele IDs that have relevant flags that have changed since the since date
        """
        return flag_ids_to(Allele, FlagComment.objects.filter(
            flag__flag_type__in={
                allele_flag_types.allele_37_not_38
            },
            created__gte=self.since
        ))

    def passes_since(self, allele_data: AlleleData) -> bool:
        """
        Is there anything about this AlleleData that indicates it should be included since the since date
        """
        if not self.since:
            return True
        if allele_data.allele_id in self.since_alleles:
            return True
        for cm in allele_data.cms:
            if cm.modified >= self.since or cm.classification.modified > self.since:
                return True
            if cm.classification_id in self.since_classifications:
                return True
        return False

    def filter_errors(self, allele_data: AlleleData) -> AlleleData:
        """
        After a since check, filter out classifications that failed validation
        TODO: store these filtered out records for CSV export with errors
        :param allele_data:
        :return:
        """
        # can already filtered out bad AlleleIDs all together
        # have to check individual classification ids
        cms = [cm for cm in allele_data.cms if not cm.classification.withdrawn and cm.classification_id not in self.bad_classification_ids]
        if bad_transcripts := self.bad_allele_transcripts.get(allele_data.allele_id):
            cms = [cm for cm in cms if cm.transcript not in bad_transcripts]

        return AlleleData(
            source=self,
            allele_id=allele_data.allele_id,
            cms=cms
        )

    @property
    def share_levels(self) -> Set[ShareLevel]:
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
        cms = ClassificationModification.objects.filter(
            is_last_published=True,
            share_level__in=self.share_levels
        )

        if not self.since:
            # can only exclude withdrawn from complete consideration
            cms = cms.exclude(classification__withdrawn=True)

        if labs := self.include_sources:
            cms.filter(classification__lab__in=labs)
        elif self.exclude_sources:
            if exclude_orgs := [org for org in self.exclude_sources if isinstance(org, Organization)]:
                cms.exclude(classification__lab__org__in=exclude_orgs)
            if exclude_labs := [lab for lab in self.exclude_sources if isinstance(lab, Lab)]:
                cms.exclude(classification__lab__in=exclude_labs)

        cms = cms.order_by('-classification__allele_id', '-classification__id')
        cms = cms.select_related('classification', 'classification__lab', 'classification__lab__organization')

        # Always safe to exclude these (unless we want them in a CSV) even with since changes
        # couldn't show these ones if we wanted to
        cms = cms.exclude(classification__allele__isnull=True).exclude(classification__variant__isnull=True)
        cms = cms.exclude(**{f'{self.c_hgvs_col}__isnull': True})

        # PERMISSION CHECK
        cms = get_objects_for_user(self.user, ClassificationModification.get_read_perm(), cms, accept_global_perms=True)
        # can't filter out transcript versions (due to 37!=38) easily here, do that later

        if self.transcript_strategy == TranscriptStrategy.REFSEQ:
            from classification.views.classification_export_view import ALISSA_ACCEPTED_TRANSCRIPTS
            acceptable_transcripts: List[Q] = [
                Q({f'classification__{self.c_hgvs_col}__startswith': tran}) for tran in ALISSA_ACCEPTED_TRANSCRIPTS
            ]
            cms = cms.filter(**{f'{self.c_hgvs_col}__startswith': 'NM'})

        return cms

    def _allele_data(self) -> Iterator[AlleleData]:
        """
        Convert the classification query set into AlleleData
        Is not filtered at this point
        """
        allele_data: Optional[AlleleData] = None
        for cm in self.cms_qs():
            allele_id = cm.classification.allele_id
            if not allele_data or allele_id != allele_data.allele_id:
                if allele_data:
                    yield allele_data
                allele_data = AlleleData(source=self, allele_id=allele_id)
            allele_data.cms.append(cm)

        if allele_data:
            yield allele_data

    def allele_data_filtered(self) -> Iterator[AlleleData]:
        """
        Main method of getting data out of ClassificationFilter
        Is only AlleleData with classifications that have passed since checks and errors
        """
        for allele_data in self._allele_data():
            if self.passes_since(allele_data):
                if filtered := self.filter_errors(allele_data):
                    yield filtered


class FileWriter:
    """
    In memory file writer, used if export create a zip file with the data potentially split across
    multiple entries
    :var file: StringIO data to write to
    :var row_count: How many lines have been written to this file
    """

    def __init__(self):
        self.file = StringIO()
        self.row_count = 0

    def __bool__(self):
        return self.row_count > 0

    def write(self, rows: List[str], count: bool = True):
        self.file.writelines(rows)
        if count:
            self.row_count += len(rows)


class ClassificationExportFormatter2(ABC):
    """
    Extend this class to export classification data into different formats
    """

    def __init__(self, filter: ClassificationFilter):
        self.filter = filter
        self.row_count = 0
        self.file_count = 0
        self.started = datetime.utcnow()

    def filename(self, part: Optional[int] = None, extension_override: Optional[str] = None) -> str:
        """
        Generate a filename
        :param part: If the data is being split, what index is this file (otherwise None)
        :param extension_override: If creating a wrapper file, e.g. "zip"
        :return: The appropriate filename
        """

        # FIXME only include genome_build if it's relevant to the file format (e.g. not JSON)
        fn = f"classifications_{self.filter.date_str}_{self.filter.genome_build}"
        if part is not None:
            # TODO make part 2 digits
            fn = f"{fn}_part_{part+1}"
        fn += f".{extension_override or self.extension()}"
        return fn

    def serve(self) -> HttpResponseBase:
        """
        Start generating the data and return it in a HTTP Response
        """
        if self.filter.row_limit:
            response = HttpResponse(content_type='application/zip')
            with zipfile.ZipFile(response, 'w') as zf:
                for index, entry in enumerate(self.yield_files()):
                    self.file_count += 1
                    # can we be more efficient than converting StringIO to a string to be converted back into bytes?
                    zf.writestr(self.filename(part=index), str(entry.file.getvalue()))
            response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
            return response
        else:
            self.file_count = 1
            # can stream in single file
            response = StreamingHttpResponse(self.yield_file(), self.content_type())
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
            return response

    @abstractmethod
    def content_type(self) -> str:
        """
        :return: Http content type
        """
        pass

    def yield_files(self) -> Iterator[FileWriter]:
        """
        :return: An iterator of chunked file data (each value could be made into a ZipEntry)
        """
        fw: Optional[FileWriter] = None

        for allele_data in self.filter.allele_data_filtered():
            to_rows = self.row(allele_data)
            self.row_count += len(to_rows)
            if not fw or (self.filter.row_limit and fw.row_count + len(to_rows) > self.filter.row_limit):
                if fw:
                    fw.write(self.footer(), count=False)
                    yield fw
                fw = FileWriter()
                fw.write(self.header(), count=False)

            fw.write(to_rows, count=True)

        if fw:
            fw.write(self.footer(), count=False)
            yield fw
        self.report_stats()

    def yield_file(self) -> Iterator[str]:
        """
        :return: An iterator for a single streaming file, call either this or yield_file
        """
        for header in self.header():
            yield header
        for allele_data in self.filter.allele_data_filtered():
            for row in self.row(allele_data):
                self.row_count += 1
                yield row
        for footer in self.footer():
            yield footer
        self.report_stats()

    @abstractmethod
    def extension(self) -> str:
        """
        Filename extension, e.g. csv, json
        """
        pass

    def header(self) -> List[str]:
        """
        Return rows to start each file, typically 0 to 1 line
        :return: A list of rows to be \n at the top of each file
        """
        return list()

    @abstractmethod
    def row(self, allele_data: AlleleData) -> List[str]:
        """
        Return the row or rows that represent this AlleleData
        :param allele_data: All the valid classifications for an Allele
        """
        pass

    def footer(self) -> List[str]:
        """
        Return rows to end each file, typically 0
        :return: A list of rows to be \n at the bottom of each file
        """
        return list()

    def report_stats(self):
        """
        TODO rename to send_report_stats
        Alerts Slack of download
        """

        # don't report bots downloading
        user = self.filter.user
        if user.groups.filter(name=bot_group().name):
            return
        end = datetime.utcnow()
        row_count = self.row_count

        body_parts = [f":simple_smile: {user.username}"]
        if request := get_current_request():
            body_parts.append(f"URL : `{request.path_info}`")
        body_parts.append(f"Rows Downloaded : *{row_count}*")
        if self.file_count > 1:
            body_parts.append(f"File Count : *{self.file_count}*")

        nb = NotificationBuilder(message="Classification Download")\
            .add_header(":arrow_down: Classification Download Completed")\
            .add_markdown("\n".join(body_parts), indented=True)
        if request := get_current_request():
            for key, value in request.GET.items():
                nb.add_field(key, value)
        nb.add_field("Duration", str((end - self.started).seconds) + " seconds")
        nb.send()