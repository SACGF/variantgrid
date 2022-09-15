from dataclasses import dataclass
from itertools import groupby
from typing import Optional, List, Iterable, Any, TypeVar, Generic, Set, Dict

from django.contrib.auth.models import User
from lazy import lazy

from classification.enums import SpecialEKeys, CriteriaEvaluation, ShareLevel
from classification.models import ClassificationModification, EvidenceKeyMap, CuratedDate, ConditionResolved, \
    classification_flag_types
from classification.models.evidence_mixin import CriteriaStrength
from classification.models.flag_types import ClassificationFlagTypes
from flags.models import Flag, FlagStatus
from genes.hgvs import CHGVS, PHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Allele, GenomeBuild, Lab

D = TypeVar("D")


@dataclass(frozen=True)
class MultiValues(Generic[D]):

    values: List[D]
    uniform: bool

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    @staticmethod
    def convert(all_value_sets: List[Set[D]]) -> 'MultiValues[D]':
        first_value_set = all_value_sets[0]
        values = list(first_value_set)
        values.sort()

        uniform = True
        for other_value_set in all_value_sets[1:]:
            if other_value_set and other_value_set != first_value_set:
                uniform = False
        return MultiValues(values=values, uniform=uniform)


@dataclass
class ClassificationGroupEntry:
    modification: ClassificationModification
    genome_build: GenomeBuild
    clinical_significance_old: Optional[str] = None
    clinical_significance_pending: Optional[str] = None

    @property
    def clin_sig(self):
        return self.modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @lazy
    def c_hgvs(self):
        return ClassificationGroup.c_hgvs_for(self.modification, self.genome_build)

    @property
    def condition_sorter(self) -> Optional[str]:
        if resolved := self.modification.classification.condition_resolution_obj:
            if resolved.terms:
                return "A" + (resolved.terms[0].name or resolved.terms[0].id).lower()
        return "Z" + (self.modification.get(SpecialEKeys.CONDITION) or "").lower()

    @lazy
    def grouping_key(self):
        clin_sig_sorter = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).classification_sorter_value
        return (
            clin_sig_sorter(self.clin_sig),
            clin_sig_sorter(self.clinical_significance_old),
            clin_sig_sorter(self.clinical_significance_pending),
            self.modification.classification.clinical_grouping_name,
            self.modification.classification.lab.organization.name,
            self.modification.classification.lab.name,
            self.c_hgvs,
            self.condition_sorter
        )

    @property
    def curated_date_check(self) -> 'CuratedDate':
        return self.modification.curated_date_check

    def __lt__(self, other):
        if self.grouping_key < other.grouping_key:
            return True


class ClassificationGroupUtils:

    def __init__(
            self,
            modifications: Optional[Iterable[ClassificationModification]] = None,
            old_modifications: Optional[Iterable[ClassificationModification]] = None,
            calculate_pending: bool = True):
        self._modifications = modifications
        self._old_modifications = old_modifications
        self.calculate_pending = calculate_pending

    @lazy
    def _pending_changes_flag_map(self) -> Dict[int, str]:
        mod_id_to_clin_sig: Dict[int, str] = dict()

        flags_qs = Flag.objects.filter(
            flag_type=classification_flag_types.classification_pending_changes,
            resolution__status=FlagStatus.OPEN)

        if self._modifications:
            flag_collections_ids = {mod.classification.flag_collection_id for mod in self._modifications}
            flags_qs.filter(collection_id__in=flag_collections_ids)

        for flag in flags_qs:
            mod_id_to_clin_sig[flag.collection_id] = flag.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY) if flag.data else 'Unknown'

        return mod_id_to_clin_sig

    @lazy
    def _classification_to_old_clin_sig(self) -> Dict[int, str]:
        old_ids: Dict[int, str] = dict()
        if self._old_modifications:
            for om in self._old_modifications:
                old_ids[om.classification_id] = om.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        return old_ids

    def pending_changes_for(self, modification: ClassificationModification) -> Optional[str]:
        if not self.calculate_pending:
            return None
        return self._pending_changes_flag_map.get(modification.classification.flag_collection_id)

    @property
    def any_pending_changes(self) -> bool:
        if not self.calculate_pending:
            return False
        return bool(self._pending_changes_flag_map)

    def map(self, modification: ClassificationModification, genome_build: GenomeBuild) -> ClassificationGroupEntry:
        clinical_significance_old = self._classification_to_old_clin_sig.get(modification.classification_id)
        clinical_significance_pending = None if not self.calculate_pending else self._pending_changes_flag_map.get(modification.classification.flag_collection_id)
        clinical_significance_current = modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        return ClassificationGroupEntry(
            modification=modification,
            genome_build=genome_build,
            clinical_significance_old=clinical_significance_old if clinical_significance_old and clinical_significance_old != clinical_significance_current else None,
            clinical_significance_pending=clinical_significance_pending if clinical_significance_pending and clinical_significance_pending != clinical_significance_current else None
        )


class ClassificationGroup:

    def __init__(self,
                 group_entries: List[ClassificationGroupEntry],
                 genome_build: GenomeBuild,
                 group_id: Optional[int] = None):

        # filter out modifications that share the same (non None) patient id
        # as in we don't want to show multiple records for the same patient
        first = group_entries[0]
        # TODO rename ClassificationGroupEntry to match
        self.clinical_significance_pending = first.clinical_significance_pending
        self.clinical_significance_old = first.clinical_significance_old
        self.group_entries = group_entries

        modification_list = [ge.modification for ge in group_entries]
        modification_list.sort(key=ClassificationGroup.sort_modifications, reverse=True)
        modification_result = list()
        seen_patient_ids = set()
        self.excluded_record_count = 0

        for modification in modification_list:
            if patient_id := modification.get(SpecialEKeys.PATIENT_ID):
                if patient_id in seen_patient_ids:
                    self.excluded_record_count += 1
                    continue
                seen_patient_ids.add(patient_id)
            modification_result.append(modification)
        self.modifications = modification_result

        self.group_id = group_id
        self.genome_build = genome_build
        self.clinical_significance_score = 0

        self.sort_order = 0  # override this once we've sorted all classifications together
        # for the sake of a JavaScript sort

    @property
    def date_sort_order_str(self) -> str:
        return str(self.sort_order).rjust(5, "0")

    @property
    def allele(self) -> Optional[Allele]:
        # method is a holdover from when allele wasn't directly accessible on a classification
        return self.most_recent.classification.allele

    def diff_ids(self) -> str:
        return ",".join([str(cm.classification_id) for cm in self.modifications])

    @property
    def gene_symbol(self) -> Optional[str]:
        return self.most_recent.get(SpecialEKeys.GENE_SYMBOL)

    @staticmethod
    def sort_modifications(mod1: ClassificationModification):
        curated_date = mod1.curated_date_check
        return [not mod1.classification.withdrawn, curated_date]

    @property
    def most_recent(self) -> ClassificationModification:
        return self.modifications[0]

    @property
    def clinical_significance(self) -> str:
        return self.most_recent.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @property
    def clinical_significance_effective(self):
        if override := self.clinical_significance_pending:
            return "pending"  # TODO, this is "pending" to help with the css class, maybe not best method name
        return self.clinical_significance or ''

    @property
    def clinical_grouping(self) -> str:
        return self.most_recent.classification.clinical_grouping_name

    @property
    def lab(self) -> Lab:
        return self.most_recent.classification.lab

    @property
    def users(self) -> List[User]:
        data = {cm.classification.user for cm in self.modifications}
        data_list = list(data)
        data_list.sort(key=lambda x: x.username)
        return data_list

    @property
    def is_withdrawn(self) -> bool:
        return all(cm.classification.withdrawn for cm in self.modifications)

    @property
    def is_discordant(self):
        if self.most_recent.share_level in ShareLevel.DISCORDANT_LEVEL_KEYS:
            if cc := self.most_recent.classification.clinical_context:
                return cc.is_discordant()

    @staticmethod
    def c_hgvs_for(cm: ClassificationModification, genome_build: GenomeBuild) -> CHGVS:
        c_parts: CHGVS
        if c_str := cm.classification.get_c_hgvs(genome_build):
            c_parts = CHGVS(c_str)
            c_parts.is_normalised = True
            c_parts.genome_build = genome_build
        else:
            c_parts = cm.classification.c_parts
            c_parts.is_normalised = False
            try:
                c_parts.genome_build = cm.classification.get_genome_build()
                c_parts.is_desired_build = genome_build.name == c_parts.genome_build.name
            except ValueError:
                pass
        return c_parts

    @lazy
    def c_hgvses(self) -> List[CHGVS]:
        unique_c = set()
        for ge in self.group_entries:
            unique_c.add(ge.c_hgvs)
        return sorted(unique_c)

    @property
    def c_hgvs(self) -> CHGVS:
        return self.c_hgvses[0]

    @lazy
    def flag_types(self) -> List[List[str]]:
        FLAG_TYPES = [
            classification_flag_types.matching_variant_flag,
            classification_flag_types.matching_variant_warning_flag,
            classification_flag_types.transcript_version_change_flag
        ]
        types: Set[str] = set()
        flags_qs = self.most_recent.classification.flag_collection_safe.flags(only_open=True)
        for type_n_res in flags_qs.filter(flag_type__in=FLAG_TYPES).values_list('flag_type', 'resolution'):
            types.add(type_n_res[0] + '.' + type_n_res[1])
        type_list = list(types)
        type_list.sort()
        return [part.split(".") for part in type_list]

    def p_hgvses(self) -> List[PHGVS]:
        unique_p = set()
        for cm in self.modifications:
            if p_hgvs := cm.p_parts:
                unique_p.add(p_hgvs.without_transcript)
        p_list = list(unique_p)
        p_list.sort()
        return p_list

    def count(self) -> int:
        return len(self.modifications)

    @lazy
    def acmg_criteria(self) -> MultiValues[CriteriaStrength]:

        def criteria_converter(cm: ClassificationModification) -> Set[CriteriaStrength]:
            strengths: Set[CriteriaStrength] = set()
            for e_key in EvidenceKeyMap.cached().criteria():
                strength = cm.get(e_key.key)
                if CriteriaEvaluation.is_met(strength):
                    strengths.add(CriteriaStrength(e_key, strength))
            return strengths

        return MultiValues.convert([criteria_converter(cm) for cm in self.modifications])

    def _evidence_key_set(self, key: str) -> List[str]:
        all_values = set()
        for cm in self.modifications:
            if value := cm.get(key):
                if isinstance(value, str):
                    value = list([value])
                all_values = all_values.union(value)
        all_values = list(all_values)
        all_values.sort()
        return all_values

    @lazy
    def zygosities(self) -> List[str]:
        return self._evidence_key_set(SpecialEKeys.ZYGOSITY)

    @lazy
    def allele_origins(self) -> List[str]:
        return self._evidence_key_set(SpecialEKeys.ALLELE_ORIGIN)

    @lazy
    def most_recent_curated(self) -> CuratedDate:
        return self.most_recent.curated_date_check

    def __lt__(self, other):
        if my_curated := self.most_recent.curated_date:
            if other_curated := other.most_recent.curated_date:
                return my_curated < other_curated
        return self.most_recent.classification.created < other.most_recent.classification.created

    def conditions(self) -> List[ConditionResolved]:
        all_terms = set()
        all_plain_texts = set()
        for cm in self.modifications:
            c = cm.classification
            if resolved := c.condition_resolution_obj:
                for term in resolved.terms:
                    all_terms.add(term)
            else:
                if text := cm.get(SpecialEKeys.CONDITION):
                    all_plain_texts.add(text)
        all_condition_resolved = list()
        for term in all_terms:
            all_condition_resolved.append(ConditionResolved(terms=[term], join=None))
        for plain_text in all_plain_texts:
            all_condition_resolved.append(ConditionResolved(terms=list(), join=None, plain_text=plain_text))

        all_condition_resolved.sort()
        return all_condition_resolved

    # def sub_groups(self) -> Optional[List['ClassificationGroup']]:
    #     if len(self.modifications) > 1:
    #         return [ClassificationGroup([cm], genome_build=self.genome_build) for cm in self.modifications]
    #     return None


class ClassificationGroups:

    def __init__(self,
                 classification_modifications: Iterable[ClassificationModification],
                 genome_build: Optional[GenomeBuild] = None,
                 group_utils: Optional[ClassificationGroupUtils] = None):

        if not genome_build:
            genome_build = GenomeBuildManager.get_current_genome_build()
        self.genome_build = genome_build

        if not group_utils:
            group_utils = ClassificationGroupUtils(classification_modifications)

        groups: List[ClassificationGroup] = list()

        classification_group_entries = [group_utils.map(m, genome_build) for m in classification_modifications]
        classification_group_entries.sort()
        for _, group in groupby(classification_group_entries, lambda x: x.grouping_key):
            group = list(group)
            actual_group = ClassificationGroup(
                group_entries=group,
                genome_build=genome_build,
                group_id=len(group) + 1,
            )
            groups.append(actual_group)
        groups.sort(key=lambda cg: cg.most_recent_curated)
        next_sort_order = 1
        for group in groups:
            group.sort_order = next_sort_order
            next_sort_order += 1
        groups.reverse()  # default store records as most recent to least recent
        self.groups = groups

    def __iter__(self):
        return iter(self.groups)

    def __len__(self):
        return len(self.groups)

    def __bool__(self):
        return bool(self.groups)

    @property
    def modifications(self) -> Iterable[ClassificationModification]:
        """
        Returns all classifications contained in this group, not just the latest
        """
        for group in self.groups:
            for record in group.modifications:
                yield record

    @property
    def latest(self) -> Iterable[ClassificationModification]:
        for group in self.groups:
            yield group.modifications[0]
