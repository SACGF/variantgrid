from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from functools import cached_property
from typing import List, Union, Iterable, Optional, Dict, Tuple, Set, Any, Mapping

from django.contrib.auth.models import User
from django.db.models import Q, Count
from django.db.models.query import QuerySet
from django.http.response import StreamingHttpResponse
from threadlocals.threadlocals import get_current_request

from classification.enums import SpecialEKeys
from classification.models.classification import ClassificationModification, \
    Classification
from classification.models.evidence_key import EvidenceKeyMap, EvidenceKey
from classification.models.flag_types import classification_flag_types
from flags.models import FlagComment
from flags.models.enums import FlagStatus
from flags.models.models import Flag
from genes.hgvs import CHGVS
from library.guardian_utils import bot_group
from library.log_utils import log_traceback, report_exc_info, report_message, NotificationBuilder
from library.utils import delimited_row, DebugTimer, local_date_string
from snpdb.models import Contig
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_genome import GenomeBuild, GenomeBuildContig
from snpdb.models.models_variant import VariantAllele, Allele


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
        self.has_explain = False


@dataclass
class KeyProperty:
    key: str
    # warning, naming this variable "property" causes confusion for the compiler below with @property
    # doesn't seem to cause any issues during runtime though
    property: str

    @property
    def field(self) -> str:
        return f"published_evidence__{self.key}__{self.property}"

    @property
    def count_key(self) -> str:
        return f"{self.field}__count"

    def count_aggregate(self) -> Count:
        return Count(f"{self.field}")

    def apply_to(self, used_key_dict: Dict[str, UsedKey]):
        used_key = used_key_dict.get(self.key)
        if not used_key:
            used_key = UsedKey()
            used_key_dict[self.key] = used_key

        if self.property == 'value':
            used_key.has_value = True
        elif self.property == 'note':
            used_key.has_note = True
        elif self.property == 'explain':
            used_key.has_explain = True


class UsedKeyTracker:

    def __init__(self,
                 user: User,
                 ekeys: EvidenceKeyMap,
                 key_value_formatter: KeyValueFormatter,
                 pretty: bool = False,
                 include_explains: bool = False,
                 ignore_evidence_keys: Optional[Set[str]] = None):
        self.user = user
        self.ekeys = ekeys
        self.key_value_formatter = key_value_formatter
        self.calc_dict: Dict[str, UsedKey] = {}
        self.pretty = pretty
        self.ordered_keys = None
        self.include_explains = include_explains
        self.ignore_evidence_keys = ignore_evidence_keys

    @property
    def keys_ignore_exclude(self) -> Iterable[EvidenceKey]:
        if self.ignore_evidence_keys:
            return [e_key for e_key in self.ekeys.all_keys if e_key.key not in self.ignore_evidence_keys]
        else:
            return self.ekeys.all_keys

    def all_key_properties(self) -> List[KeyProperty]:
        all_props: List[KeyProperty] = []
        properties = ['value', 'note']
        if self.include_explains:
            properties.append('explain')
        for e_key in self.keys_ignore_exclude:
            for prop in properties:
                all_props.append(KeyProperty(key=e_key.key, property=prop))
        return all_props

    def check_evidence_qs(self, qs: QuerySet[ClassificationModification]):
        all_key_properties = self.all_key_properties()
        aggregate_list = [kp.count_aggregate() for kp in all_key_properties]
        result_dict = qs.aggregate(*aggregate_list)
        for kp in all_key_properties:
            if result_dict.get(kp.count_key):
                kp.apply_to(self.calc_dict)

    def check_record(self, vcm: ClassificationModification):
        self.check_evidence(vcm.evidence)

    def check_evidence(self, evidence: Dict[str, Any]):
        has_value = False
        has_note = False
        has_explain = False
        for key, valueObj in evidence.items():
            if isinstance(valueObj, Mapping):
                has_value = valueObj.get('value') is not None
                has_note = valueObj.get('note') is not None
                has_explain = valueObj.get('explain') is not None

            if has_value or has_note or (self.include_explains and has_explain):
                used_key = self.calc_dict.get(key)
                if not used_key:
                    used_key = UsedKey()
                    self.calc_dict[key] = used_key

                used_key.has_value = used_key.has_value or has_value
                used_key.has_note = used_key.has_note or has_note
                used_key.has_explain = used_key.has_explain or has_explain

    def process(self):
        self.ordered_keys = []
        for ekey in self.ekeys.all_keys:
            used_key = self.calc_dict.get(ekey.key)
            if used_key:
                used_key.ekey = ekey
                self.ordered_keys.append(used_key)

    def header(self) -> List[str]:
        self.process()
        cols: List[str] = []
        for used_key in self.ordered_keys:
            if used_key.has_value:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, pretty=self.pretty))
            if used_key.has_note:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, is_note=True, pretty=self.pretty))
            if self.include_explains and used_key.has_explain:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, pretty=self.pretty) + '.explain')
        return cols

    def row(self, classification_modification: ClassificationModification) -> List[Optional[str]]:
        cols: List[Optional[str]] = []
        evidence = classification_modification.get_visible_evidence(self.user)
        for used_key in self.ordered_keys:
            value_obj = evidence.get(used_key.ekey.key)
            if used_key.has_value:
                if not value_obj:
                    cols.append(None)
                else:
                    value = value_obj.get('value')
                    cols.append(self.key_value_formatter.value_for(used_key.ekey, value, pretty=self.pretty))
            if used_key.has_note:
                if not value_obj:
                    cols.append(None)
                else:
                    cols.append(value_obj.get('note'))
            if self.include_explains and used_key.has_explain:
                if not value_obj:
                    cols.append(None)
                else:
                    cols.append(value_obj.get('explain'))

        return cols


class ConflictStrategy(str, Enum):
    MOST_BENIGN = 'most_benign'
    MOST_PATHOGENIC = 'most_pathogenic'


class VCFEncoding:
    BASIC = 'basic'
    FULL = 'full'


class VariantWithChgvs:

    def __init__(self, vcm: ClassificationModification, chgvs: CHGVS):
        self.vcm = vcm
        self.chgvs = chgvs

    @cached_property
    def transcript_version(self) -> int:
        if self.chgvs.transcript_parts:
            return self.chgvs.transcript_parts or 0
        else:
            return 0

    @property
    def c_hgvs_without_transcript_version(self) -> CHGVS:
        return self.chgvs.without_transcript_version


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
    def different_c_hgvs(self):
        first_c = self.vcmcs[0].chgvs.without_transcript_version
        for vc in self.vcmcs[1:]:
            if first_c != vc.chgvs.without_transcript_version:
                return True
        return False

    @property
    def cms(self) -> List[ClassificationModification]:
        return [vcmcs.vcm for vcmcs in self.vcmcs]

    @property
    def chgvs(self):
        return self.highest_transcript_chgvs


class AlleleGroup:
    """
    A bunch of records linked to variants with the same allele
    Also contains a "target_variant" from the desired genome_build
    """

    def __init__(self, source: 'ExportFormatter', allele_id: int, allele_flag_collection_id: int, genome_build: GenomeBuild):
        self.source = source
        self.allele_id = allele_id
        self.allele_flag_collection_id = allele_flag_collection_id
        self.target_variant = None
        self.genome_build = genome_build
        self.variant_ids = []
        self.data: List[ClassificationModification] = []
        self.withdrawn: List[ClassificationModification] = []
        self.failed: List[ClassificationModification] = []
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

    def iter_c_hgvs_versionless_transcripts(self) -> Iterable[Tuple[CHGVS, List[VariantWithChgvs]]]:
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
