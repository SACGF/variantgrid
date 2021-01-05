import collections
from collections import defaultdict
from datetime import datetime

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models.query import QuerySet
from django.http.response import StreamingHttpResponse
from typing import List, Union, Iterable, Optional, Dict, Tuple, Set

from django.utils.timezone import now
from lazy import lazy
from pytz import timezone

from flags.models import FlagComment
from flags.models.enums import FlagStatus
from flags.models.models import Flag
from genes.hgvs import CHGVS
from library.log_utils import log_traceback, report_exc_info, report_message
from library.utils import delimited_row
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_genome import GenomeBuild, GenomeBuildContig
from snpdb.models.models_variant import VariantAllele, Allele
from classification.enums import SpecialEKeys
from classification.models.evidence_key import EvidenceKeyMap, EvidenceKey
from classification.models.flag_types import classification_flag_types
from classification.models.classification import ClassificationModification,\
    Classification


class KeyValueFormatter:

    def header_for(self, ekey: EvidenceKey, is_note: bool = False, pretty: bool = False) -> str:
        label: str
        if pretty:
            label = ekey.pretty_label
            if is_note:
                label += ' note'
        else:
            label = ekey.key
            if is_note:
                label += '.note'
        return label

    def value_for(self, ekey: EvidenceKey, value, pretty: bool = False):
        if pretty:
            return ekey.pretty_value(value)

        if isinstance(value, list):
            value = ', '.join((str(item) for item in value))
        elif value is True:
            return 'TRUE'
        elif value is False:
            return 'FALSE'
        return value


class UsedKey:

    def __init__(self):
        self.ekey = None
        self.has_value = False
        self.has_note = False


class UsedKeyTracker:

    def __init__(self, user: User, ekeys: EvidenceKeyMap, key_value_formatter: KeyValueFormatter, pretty: bool = False):
        self.user = user
        self.ekeys = ekeys
        self.key_value_formatter = key_value_formatter
        self.calc_dict = {}
        self.pretty = pretty
        self.ordered_keys = None

    def check_record(self, vcm: ClassificationModification):
        has_value = False
        has_note = False
        for key, valueObj in vcm.evidence.items():
            if isinstance(valueObj, collections.Mapping):
                has_value = valueObj.get('value') is not None
                has_note = valueObj.get('note') is not None

            if has_value or has_note:
                used_key = self.calc_dict.get(key)
                if not used_key:
                    used_key = UsedKey()
                    self.calc_dict[key] = used_key

                used_key.has_value = used_key.has_value or has_value
                used_key.has_note = used_key.has_note or has_note

    def process(self):
        self.ordered_keys = []
        for ekey in self.ekeys.all_keys:
            used_key = self.calc_dict.get(ekey.key)
            if used_key:
                used_key.ekey = ekey
                self.ordered_keys.append(used_key)

    def header(self) -> List[str]:
        self.process()
        cols = []
        for used_key in self.ordered_keys:
            if used_key.has_value:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, pretty=self.pretty))
            if used_key.has_note:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, is_note=True, pretty=self.pretty))
        return cols

    def row(self, classification_modification: ClassificationModification) -> List:
        cols = []
        evidence = classification_modification.get_visible_evidence(self.user)
        for used_key in self.ordered_keys:
            valueObj = evidence.get(used_key.ekey.key)
            if used_key.has_value:
                if not valueObj:
                    cols.append(None)
                else:
                    value = valueObj.get('value')
                    cols.append(self.key_value_formatter.value_for(used_key.ekey, value, pretty=self.pretty))
            if used_key.has_note:
                if not valueObj:
                    cols.append(None)
                else:
                    cols.append(valueObj.get('note'))
        return cols


class ConflictStrategy:
    MOST_BENIGN = 'most_benign'
    MOST_PATHOGENIC = 'most_pathogenic'

class VCFEncoding:
    BASIC = 'basic'
    FULL = 'full'

class VariantWithChgvs:

    def __init__(self, vcm: ClassificationModification, chgvs: CHGVS):
        self.vcm = vcm
        self.chgvs = chgvs

    @lazy
    def transcript_version(self) -> int:
        if self.chgvs.transcript_parts:
            return self.chgvs.transcript_parts or 0
        else:
            return 0

    @lazy
    def c_hgvs_without_transcript_version(self) -> CHGVS:
        return self.chgvs.without_transcript_version()


class TranscriptGroup:

    def __init__(self):
        self.highest_transcript_version: Optional[int] = None
        self.highest_transcript_chgvs: Optional[CHGVS] = None
        self.vcmcs: List[VariantWithChgvs] = []

    def add(self, vcmc: VariantWithChgvs):
        self.vcmcs.append(vcmc)

        transcript_version = vcmc.transcript_version
        if self.highest_transcript_version is None or transcript_version > self.highest_transcript_version:
            self.highest_transcript_version = transcript_version
            self.highest_transcript_chgvs = vcmc.chgvs

    @property
    def chgvs(self):
        return self.highest_transcript_chgvs


class AlleleGroup:
    """
    A bunch of records linked to variants with the same allele
    Also contains a "target_variant" from the desired genome_build
    """

    def __init__(self, source: 'ExportFormatter', allele_id: int, genome_build: GenomeBuild):
        self.source = source
        self.allele_id = allele_id
        self.target_variant = None
        self.genome_build = genome_build
        self.variant_ids = []
        self.data: List[ClassificationModification] = []
        self.source = source

    def filter_out_transcripts(self, transcripts: Set[str]) -> List[ClassificationModification]:
        """
        Returns ClassificationModificats there weren't included due to errors
        """
        passes: List[ClassificationModification] = []
        fails: List[ClassificationModification] = []

        for vcm in self.data:
            if vcm.transcript in transcripts:
                fails.append(vcm)
            else:
                passes.append(vcm)
        if fails:
            self.data = passes

        return fails

    def liftover(self, vcm: ClassificationModification) -> CHGVS:
        chgvs_str = vcm.classification.get_c_hgvs(self.genome_build, use_full=self.source.use_full_chgvs)
        if not chgvs_str:
            raise ValueError(f"Unable to generate c.hgvs full={self.source.use_full_chgvs}, this record should have been filtered out")
        return CHGVS(vcm.classification.get_c_hgvs(self.genome_build, use_full=self.source.use_full_chgvs))

    def iter_c_hgvs(self) -> Iterable[Tuple[CHGVS, List[ClassificationModification]]]:
        by_transcript: Dict[CHGVS, List[ClassificationModification]] = defaultdict(list)

        for vcm in self.data:
            c_parts = self.liftover(vcm)
            if c_parts:
                by_transcript[c_parts].append(vcm)

        for c_hgvs, vcms in by_transcript.items():
            yield c_hgvs, vcms

    def iter_c_hgvs_versionless_transcripts(self) -> Iterable[Tuple[CHGVS, VariantWithChgvs]]:
        by_versionless_transcript: Dict[str, TranscriptGroup] = defaultdict(TranscriptGroup)

        for vcm in self.data:
            c_parts = self.liftover(vcm)
            if c_parts:
                transcript_parts = c_parts.transcript_parts
                if transcript_parts:
                    transcript_no_version = transcript_parts.identifier
                    by_versionless_transcript[transcript_no_version].add(VariantWithChgvs(vcm=vcm, chgvs=c_parts))
                else:
                    report_message('MVL export : Could not extract transcript from c.hgvs', extra_data={'chgvs': c_parts.full_c_hgvs})
            else:
                report_message('MVL export : Could not liftover', extra_data={'imported_chgvs': vcm.get(SpecialEKeys.C_HGVS), 'id': vcm.classification_id})

        for _, transcript_groups in by_versionless_transcript.items():
            yield transcript_groups.highest_transcript_chgvs, transcript_groups.vcmcs


class BaseExportFormatter:

    def __init__(self):
        pass

    def export(self, as_attachment: bool = True) -> StreamingHttpResponse:
        raise NotImplementedError("export_raw has not been implemented")


class ExportFormatter(BaseExportFormatter):
    """
    Given a stream of AlleleGroups, output string data for use in exports
    """
    @property
    def version(self):
        return '1.0'

    def __init__(self, genome_build: GenomeBuild, qs: QuerySet, user: User = None, since: datetime = None):
        self.genome_build = genome_build
        self.used_contigs = set()
        self.user = user
        self.allele_groups = []
        self.since = since

        self.error_message_ids = dict()
        self.raw_qs = qs
        self.started = datetime.utcnow()
        self.qs = qs.filter(**{f'{self.preferred_chgvs_column}__isnull': False})

        super().__init__()

    @property
    def use_full_chgvs(self) -> bool:
        """ override to return True if you want dels, ins etc to list all nucleotides after """
        return False

    @lazy
    def preferred_chgvs_column(self) -> str:
        return ClassificationModification.column_name_for_build(self.genome_build, use_full=self.use_full_chgvs)

    @lazy
    def ekeys(self) -> EvidenceKeyMap:
        return EvidenceKeyMap.cached()

    def generate_filename(self,
        prefix: str = 'classifications',
        include_date=True,
        include_genome_build: bool = True,
        suffix: str = None,
        extension: str = 'csv') -> str:

        parts = []
        if prefix:
            parts.append(prefix)
        if include_date:
            parts.append( now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d") )
        if include_genome_build:
            parts.append(self.genome_build.name)
        if suffix:
            parts.append(suffix)

        return '_'.join(parts) + '.' + extension

    def record_and_error(self):
        for vcm in ClassificationModification.objects.filter(pk__in=self.error_message_ids.keys()).order_by('created'):
            message = self.error_message_ids[vcm.id]
            yield vcm, message

    def filter_mismatched_transcripts(self, allele_group: AlleleGroup):
        flag_collection_id = Allele.objects.filter(pk=allele_group.allele_id).values_list('flag_collection_id',
                                                                                          flat=True).first()
        transcripts: Set[str] = set()
        if flag_collection_id:
            outstanding_mismatches = Flag.objects.filter(
                collection_id=flag_collection_id,
                flag_type=allele_flag_types.allele_37_not_38,
                resolution__status=FlagStatus.OPEN
            )
            outstanding_mismatch: Flag
            for outstanding_mismatch in outstanding_mismatches:
                transcript = (outstanding_mismatch.data or {}).get('transcript')
                if transcript:
                    transcripts.add(transcript)

        if transcripts:
            failed_ids = [vcm.id for vcm in allele_group.filter_out_transcripts(transcripts)]
            self.record_errors(failed_ids, f'transcript across genome builds requires confirmation')

    def passes_since_check(self, ag: AlleleGroup) -> bool:
        if self.since:
            cm: ClassificationModification
            for cm in ag.data:
                if cm.modified > self.since:
                    return True

            # no modifications passed the check, but maybe a flag has changed on the allele of classifications
            flag_collection_ids = list()
            flag_collection_ids.append( Allele.objects.values_list("flag_collection_id", flat=True).first() )
            for cm in ag.data:
                flag_collection_ids.append(cm.classification.flag_collection_id)

            if FlagComment.objects.filter(flag__collection_id__in=flag_collection_ids, created__gt=self.since).exists():
                return True

        else:
            return True


    def passes_flag_check(self, vcm: ClassificationModification) -> bool:
        outstanding_warning = Flag.objects.filter(
            collection_id=vcm.classification.flag_collection_id,
            flag_type__in=[classification_flag_types.matching_variant_warning_flag, classification_flag_types.transcript_version_change_flag],
            resolution__status=FlagStatus.OPEN
        ).exists()
        if outstanding_warning:
            self.record_errors(vcm.id, 'Requires confirmation of variant match')
            return False

        return True

    def record_errors(self, ids: Union[QuerySet, int, List[int], ClassificationModification], message: str):
        if isinstance(ids, QuerySet):
            ids = ids.values_list('id', flat=True)
        elif isinstance(ids, int):
            ids = [ids]
        elif isinstance(ids, ClassificationModification):
            ids = [ids.id]

        for record_id in ids:
            self.error_message_ids[record_id] = message

    @staticmethod
    def write_single_row(data: list, delimiter: str = ',') -> str:
        return delimited_row(data, delimiter)

    def prepare_groups(self):
        """
        Given a QuerySet, works out what variants and alleles will be involved, and then returns AlleleGroups
        containing all the relevant rows from qs, as well as the relevant variant from genome_build.
        The order of the AlleleGroups will be in locus position order of the "target_variant" (variant of the genome build),
        if no target_variant could be found for the AlleleGroup it is omitted.
        """
        could_not_liftover_qs = self.raw_qs.filter(**{f'{self.preferred_chgvs_column}__isnull': True})
        self.record_errors(could_not_liftover_qs, 'Could not liftover/normalise')

        all_variants = self.qs.values_list('classification__variant', flat=True).distinct()
        all_alleles = VariantAllele.objects.filter(variant__in=all_variants).values_list('allele', flat=True)
        all_variant_alleles = VariantAllele.objects.filter(allele__in=all_alleles).order_by('allele_id')\
            .select_related('genome_build', 'variant', 'variant__locus', 'variant__locus__contig')

        allele_group = None
        allele_warning_variant_ids = set()

        def process_allele_group(group: AlleleGroup):
            if group:
                if not group.target_variant:
                    report_message(f'Could not find "{self.genome_build}" variant for allele {group.allele_id}')
                else:
                    self.allele_groups.append(group)

        for va in all_variant_alleles:
            # Allele has changed, start a new group
            if not allele_group or allele_group.allele_id != va.allele_id:
                process_allele_group(allele_group)
                allele_group = AlleleGroup(source=self, allele_id=va.allele_id, genome_build=self.genome_build)

            if va.genome_build.is_equivalent(self.genome_build):
                variant = va.variant
                allele_group.target_variant = variant
                self.used_contigs.add(variant.locus.contig)

            allele_group.variant_ids.append(va.variant_id)

        process_allele_group(allele_group)

        contig_order: Dict[int, int] = dict()
        for cbc in GenomeBuildContig.objects.filter(genome_build=self.genome_build).order_by('order'):
            contig = cbc.contig
            contig_order[contig.id] = len(contig_order)

        def variant_position_sorter(ag: AlleleGroup):
            locus = ag.target_variant.locus
            return contig_order.get(locus.contig_id, 0), locus.position

        self.allele_groups.sort(key=variant_position_sorter)

        if allele_warning_variant_ids:
            allele_warning_variant_qs = self.qs.filter(classification__variant__in=allele_warning_variant_ids)
            self.record_errors(allele_warning_variant_qs, f'Allele has an outstanding warning')

    def iter_group_by_allele(self):
        for allele_group in self.allele_groups:
            # actually populate the allele group data now
            allele_group.data = [vcm for vcm in self.qs.filter(classification__variant__in=allele_group.variant_ids) if self.passes_flag_check(vcm)]
            self.filter_mismatched_transcripts(allele_group)

            if allele_group.data and self.passes_since_check(allele_group):
                yield allele_group

    def row_iterator(self) -> Iterable[AlleleGroup]:
        return self.iter_group_by_allele()

    def header(self) -> Optional[str]:
        """
        Header at the start of the output, (must provide your own new line characters).
        Return None for no header
        """
        return None

    def row(self, group) -> Optional[str]:
        """
        Return data for an AlleleGroup, you can return multiple lines, a single line or None
        """
        return None

    def footer(self) -> Optional[str]:
        """
        Return any row to go after all other rows
        """
        return None

    def content_type(self) -> str:
        return 'text/csv'

    def filename(self) -> str:
        return 'classifications.txt'

    def is_discordant(self, vc: Classification):
        if vc.flag_collection_id:
            return Flag.objects.filter(
                collection=vc.flag_collection_id,
                flag_type=classification_flag_types.discordant,
                resolution__status=FlagStatus.OPEN
            ).exists()
        return False

    def export(self, as_attachment: bool = True) -> StreamingHttpResponse:

        def iter_allele_to_row():
            self.prepare_groups()

            header = self.header()
            if header:
                yield header
            for vdata in self.row_iterator():
                try:
                    row = self.row(vdata)
                    if row:
                        yield row
                except GeneratorExit:
                    # user has cancelled the download, just stop now
                    return
                except:
                    print('Excepting during export')
                    log_traceback()
                    report_exc_info()

            footer = self.footer()
            if footer:
                yield footer

        response = StreamingHttpResponse(iter_allele_to_row(), content_type=self.content_type())
        modified_str = self.started.strftime("%a, %d %b %Y %H:%M:%S GMT")  # e.g. 'Wed, 21 Oct 2015 07:28:00 GMT'

        response['Last-Modified'] = modified_str
        if as_attachment:
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
        return response
