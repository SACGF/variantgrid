from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Iterable, Optional, Any, Mapping, Callable

from django.contrib.auth.models import User
from django.db.models import Count
from django.db.models.query import QuerySet

from classification.enums import CriteriaEvaluation
from classification.models.classification import ClassificationModification
from classification.models.evidence_key import EvidenceKeyMap, EvidenceKey
from genes.hgvs import CHGVS
from library.cache import clear_cached_property


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

    def value_for(self, ekey: EvidenceKey, value, pretty: bool = False, cell_formatter: Optional[Callable[[Any], Any]] = False):
        if pretty:
            value = ekey.pretty_value(value)
        else:
            if isinstance(value, list):
                value = ', '.join((str(item) for item in value))
            elif value is True:
                return 'TRUE'
            elif value is False:
                return 'FALSE'

        if cell_formatter:
            value = cell_formatter(value)

        return value


class UsedKey:

    def __init__(self):
        self.ekey: Optional[EvidenceKey] = None
        self.has_value = False
        self.has_note = False
        self.has_explain = False


@dataclass
class KeyProperty:
    key: str
    # warning, naming this variable "property" causes confusion for the compiler below with @property
    # doesn't seem to cause any issues during runtime though
    prop: str

    @property
    def id(self):
        return f"{self.key}.{self.prop}"

    @property
    def field(self) -> str:
        return f"published_evidence__{self.key}__{self.prop}"

    @property
    def count_key(self) -> str:
        return f"{self.field}__count"

    def count_aggregate(self) -> Count:
        return Count(f"{self.field}")

    def apply_to(self, used_key_dict: dict[str, UsedKey]):
        used_key = used_key_dict.get(self.key)
        if not used_key:
            used_key = UsedKey()
            used_key_dict[self.key] = used_key

        if self.prop == 'value':
            used_key.has_value = True
        elif self.prop == 'note':
            used_key.has_note = True
        elif self.prop == 'explain':
            used_key.has_explain = True


class UsedKeyTracker:

    def __init__(self,
                 user: User,
                 ekeys: EvidenceKeyMap,
                 key_value_formatter: KeyValueFormatter,
                 pretty: bool = False,
                 cell_formatter: Optional[Callable[[Any], Any]] = None,
                 include_explains_and_notes: bool = False,
                 ignore_evidence_keys: Optional[set[str]] = None,
                 include_only_evidence_keys: Optional[set[str]] = None):
        self.user = user
        self.ekeys = ekeys
        self.key_value_formatter = key_value_formatter
        self.calc_dict: dict[str, UsedKey] = {}
        self.pretty = pretty
        self.cell_formatter = cell_formatter
        self.include_explains_and_notes = include_explains_and_notes
        self.ignore_evidence_keys = ignore_evidence_keys
        self.include_only_evidence_keys = include_only_evidence_keys

    @property
    def considered_keys(self) -> Iterable[EvidenceKey]:
        consider_keys = self.ekeys.all_keys
        if include_only_keys := self.include_only_evidence_keys:
            consider_keys = [e_key for e_key in consider_keys if e_key.key in include_only_keys]

        if self.ignore_evidence_keys:
            consider_keys = [e_key for e_key in consider_keys if e_key.key not in self.ignore_evidence_keys]

        return consider_keys

    @cached_property
    def all_key_properties(self) -> list[KeyProperty]:
        all_props: list[KeyProperty] = []
        properties = ['value']
        if self.include_explains_and_notes:
            properties += ['note', 'explain']
        for e_key in self.considered_keys:
            for prop in properties:
                all_props.append(KeyProperty(key=e_key.key, prop=prop))
        return all_props

    @cached_property
    def key_property_map(self) -> dict[str, KeyProperty]:
        key_prop_map: dict[str, KeyProperty] = {}
        for key in self.all_key_properties:
            key_prop_map[key.id] = key
        return key_prop_map

    def check_evidence_qs(self, qs: QuerySet[ClassificationModification]):
        """
        By performing an aggregate over the entire queryset, work out which columns
        have values, the old way was to just review all the data twice.
        """
        aggregate_list = [kp.count_aggregate() for kp in self.all_key_properties]
        result_dict = qs.aggregate(*aggregate_list)
        for kp in self.all_key_properties:
            if result_dict.get(kp.count_key):
                kp.apply_to(self.calc_dict)

    def check_evidence_enable_all_considered(self):
        for kp in self.all_key_properties:
            kp.apply_to(self.calc_dict)

    def check_record(self, vcm: ClassificationModification):
        self.check_evidence(vcm.evidence)

    def check_evidence(self, evidence: dict[str, Any]):
        clear_cached_property(self, "ordered_keys")
        for key, valueObj in evidence.items():
            if isinstance(valueObj, Mapping):
                for prop in ["value", "note", "explain"]:
                    if valueObj.get(prop) is not None:
                        if key_prop := self.key_property_map.get(f"{key}.{prop}"):
                            key_prop.apply_to(self.calc_dict)

    @cached_property
    def ordered_keys(self) -> list[UsedKey]:
        ordered_keys: list[UsedKey] = []
        for ekey in self.considered_keys:
            used_key = self.calc_dict.get(ekey.key)
            if used_key:
                used_key.ekey = ekey
                ordered_keys.append(used_key)
        return ordered_keys

    def header(self) -> list[str]:
        cols: list[str] = []
        for used_key in self.ordered_keys:
            if used_key.has_value:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, pretty=self.pretty))
            if used_key.has_note:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, is_note=True, pretty=self.pretty))
            if used_key.has_explain:
                cols.append(self.key_value_formatter.header_for(used_key.ekey, pretty=self.pretty) + '.explain')
        return cols

    def row(self, classification_modification: ClassificationModification, formatter: Callable[[Any], Any]) -> list[Optional[str]]:
        cols: list[Optional[str]] = []
        evidence = classification_modification.get_visible_evidence(self.user)
        for used_key in self.ordered_keys:
            value_obj = evidence.get(used_key.ekey.key)
            if used_key.has_value:
                if not value_obj:
                    cols.append(None)
                else:
                    value = value_obj.get('value')

                    # prettyValue for a horak code will say Met: 2 Points
                    # we just want to say 2
                    if used_key.ekey.crit_uses_points:
                        points = CriteriaEvaluation.POINTS.get(value)
                        if points is not None:
                            value = str(points)

                    cols.append(self.key_value_formatter.value_for(used_key.ekey, value, pretty=self.pretty, cell_formatter=self.cell_formatter))

            # an earlier check determines if we even think about adding has_note or has_explain
            if used_key.has_note:
                if not value_obj:
                    cols.append(None)
                else:
                    cols.append(value_obj.get('note'))
            if used_key.has_explain:
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
        if tp := self.chgvs.transcript_parts:
            return tp.version or 0
        else:
            return 0

    @property
    def c_hgvs_without_transcript_version(self) -> CHGVS:
        return self.chgvs.without_transcript_version


class TranscriptGroup:

    def __init__(self):
        self.highest_transcript_version: Optional[int] = None
        self.highest_transcript_chgvs: Optional[CHGVS] = None
        self.vcmcs: list[VariantWithChgvs] = []

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
    def cms(self) -> list[ClassificationModification]:
        return [vcmcs.vcm for vcmcs in self.vcmcs]

    @property
    def chgvs(self):
        return self.highest_transcript_chgvs
