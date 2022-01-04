import zipfile
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from io import StringIO
from typing import Optional, Set, Union, List, Iterator, Dict, Tuple, Type, Iterable
from django.conf import settings
from django.db.models import QuerySet, Q
from django.utils.timezone import now
from pytz import timezone
from cyvcf2.cyvcf2 import defaultdict
from django.http import HttpRequest, HttpResponse, StreamingHttpResponse
from django_messages.admin import User
from guardian.shortcuts import get_objects_for_user
from lazy import lazy

from classification.enums import ShareLevel
from classification.models import ClassificationModification, classification_flag_types, Classification
from flags.models import Flag, FlagStatus, FlagComment, FlagsMixin
from genes.hgvs import CHGVS
from snpdb.models import Lab, Organization, GenomeBuild, allele_flag_types, Allele


@dataclass
class AlleleData:
    source: 'ClassificationFilter'
    allele_id: int
    cms: List[ClassificationModification] = field(default_factory=list)

    def __bool__(self):
        return bool(self.cms)


@dataclass
class CHGVSData:
    source: AlleleData
    chgvs: CHGVS
    different_versions: bool = False
    cms: List[ClassificationModification] = field(default_factory=list)

    @staticmethod
    def split_into_c_hgvs(
            allele_data: AlleleData,
            use_full: bool) -> List['CHGVSData']:
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


def flag_ids_to(model: Type[FlagsMixin], qs: QuerySet) -> Set[int]:
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
    ALL = "all"
    REFSEQ = "refseq"


@dataclass
class ClassificationFilter:
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
        return now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")

    @staticmethod
    def string_to_group_name(model: Type[Union[Lab, Organization]], group_names: str) -> Set:
        if not group_names:
            return set()
        parts = [gn.strip() for gn in group_names.split(',')]
        objs = [model.objects.filter(group_name=group_name).first() for group_name in parts]
        objs = [m for m in objs if m]
        return set(objs)

    @staticmethod
    def from_request(request: HttpRequest) -> 'ClassificationFilter':
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

        # FIXME include row_limit into filter

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
        if self.genome_build == GenomeBuild.grch37():
            return 'classification__chgvs_grch37'
        else:
            return 'classification__chgvs_grch38'

    @lazy
    def bad_allele_transcripts(self) -> Dict[int, Set[str]]:
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
        return flag_ids_to(Classification, Flag.objects.filter(
            flag_type__in={classification_flag_types.transcript_version_change_flag,
                           classification_flag_types.matching_variant_warning_flag},
            resolution__status=FlagStatus.OPEN
        ))

    @lazy
    def discordant_classification_ids(self) -> Set[int]:
        return flag_ids_to(Classification, Flag.objects.filter(
            flag_type=classification_flag_types.discordant,
            resolution__status=FlagStatus.OPEN
        ))

    def is_discordant(self, cm: ClassificationModification):
        return cm.classification_id in self.discordant_classification_ids

    @lazy
    def since_classifications(self) -> Set[int]:
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
        return flag_ids_to(Allele, FlagComment.objects.filter(
            flag__flag_type__in={
                allele_flag_types.allele_37_not_38
            },
            created__gte=self.since
        ))

    def passes_since(self, allele_data: AlleleData):
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

    def filter_errors(self, allele_data: AlleleData):
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
    def share_levels(self):
        share_levels: Set[ShareLevel] = set()
        for sl in ShareLevel.ALL_LEVELS:
            if sl >= self.min_share_level:
                share_levels.add(sl)
        return share_levels

    def _allele_data(self) -> Iterator[AlleleData]:
        """
        Given all the parameters, as efficiently as it can,
        :return:
        """
        print("SHARE LEVELS = ")
        print(self.share_levels)
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

        allele_data: Optional[AlleleData] = None
        for cm in cms:
            allele_id = cm.classification.allele_id
            if not allele_data or allele_id != allele_data.allele_id:
                if allele_data:
                    yield allele_data
                allele_data = AlleleData(source=self, allele_id=allele_id)
            allele_data.cms.append(cm)

        if allele_data:
            yield allele_data

    def allele_data_filtered(self):
        for allele_data in self._allele_data():
            if self.passes_since(allele_data):
                if filtered := self.filter_errors(allele_data):
                    yield filtered


class FileWriter:

    def __init__(self):
        self.file = StringIO()
        self.row_count = 0

    def __bool__(self):
        return self.row_count > 0

    def write(self, rows: List[str], count: bool = True):
        self.file.writelines(rows)
        if count:
            self.row_count += len(rows)


class ClassificationExportFormatter2:

    def __init__(self, filter: ClassificationFilter):
        self.filter = filter

    def filename(self, part: Optional[int] = None, extension_override: Optional[str] = None):
        # FIXME only include genome_build if it's relevant to the file format (e.g. not JSON)
        fn = f"classifications_{self.filter.date_str}_{self.filter.genome_build}"
        if part is not None:
            # TODO make part 2 digits
            fn = f"{fn}_part_{part+1}"
        fn += f".{extension_override or self.extension()}"
        return fn

    def serve(self):
        if self.filter.row_limit:
            response = HttpResponse(content_type='application/zip')
            with zipfile.ZipFile(response, 'w') as zf:
                for index, entry in enumerate(self.yield_files()):
                    # can we be more efficient than converting StringIO to a string to be converted back into bytes?
                    zf.writestr(self.filename(part=index), str(entry.file.getvalue()))
            response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
            return response
        else:
            # can stream in single file
            response = StreamingHttpResponse(self.yield_file(), self.content_type())
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
            return response

    def content_type(self) -> str:
        return "text/csv"

    def yield_files(self) -> Iterator[FileWriter]:
        """
        :return: An iterator of chunked file data (each value could be made into a ZipEntry)
        """
        fw: Optional[FileWriter] = None

        for allele_data in self.filter.allele_data_filtered():
            to_rows = self.row(allele_data)
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

    def yield_file(self) -> Iterator[str]:
        """
        :return: An iterator for a single streaming file
        """
        for header in self.header():
            yield header
        for allele_data in self.filter.allele_data_filtered():
            for row in self.row(allele_data):
                yield row
        for footer in self.footer():
            yield footer

    def extension(self) -> str:
        raise NotImplementedError("extension not implemented")

    def header(self) -> List[str]:
        return list()

    def row(self, allele_data: AlleleData) -> List[str]:
        raise NotImplementedError("row not implemented")

    def footer(self) -> List[str]:
        return list()
