from dataclasses import dataclass
from itertools import groupby
from typing import Optional, List, Iterable, Any, TypeVar, Generic, Set, Dict

from django.contrib.auth.models import User
from lazy import lazy

from classification.enums import SpecialEKeys, CriteriaEvaluation, ShareLevel
from classification.models import ClassificationModification, EvidenceKeyMap, CuratedDate, ConditionResolved, \
    classification_flag_types
from classification.models.evidence_mixin import CriteriaStrength
from genes.hgvs import CHGVS, PHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Allele, GenomeBuild

D = TypeVar("D")


@dataclass(frozen=True)
class MultiValues(Generic[D]):

    values: Iterable[D]
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


class ClassificationGroup:

    def __init__(self,
                 modifications: Iterable[ClassificationModification],
                 genome_build: GenomeBuild,
                 group_id: Optional[int] = None):

        # filter out modifications that share the same (non None) patient id
        # as in we don't want to show multiple records for the same patient
        modification_list = list(modifications)
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
    def sort_order_str(self) -> str:
        return str(self.sort_order).rjust(5, "0")

    @property
    def allele(self) -> Optional[Allele]:
        if variant := self.most_recent.classification.variant:
            return variant.allele
        return None

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
    def clinical_grouping(self) -> str:
        return self.most_recent.classification.clinical_grouping_name

    @property
    def organization(self) -> str:
        return self.most_recent.classification.lab.organization.name

    @property
    def lab(self) -> str:
        return self.most_recent.classification.lab.name

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
            except KeyError:
                pass
        return c_parts

    @lazy
    def c_hgvses(self) -> List[CHGVS]:
        unique_c = set()
        for cm in self.modifications:
            unique_c.add(ClassificationGroup.c_hgvs_for(cm, self.genome_build))
        c_list = list(unique_c)
        c_list.sort()
        return c_list

    @property
    def c_hgvs(self) -> CHGVS:
        return self.c_hgvses[0]

    @lazy
    def flag_types(self) -> Set[str]:
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

        def criteria_converter(cm: ClassificationModification) -> Set[str]:
            strengths: Set[str] = set()
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

    def sub_groups(self) -> Optional[List['ClassificationGroup']]:
        if len(self.modifications) > 1:
            return [ClassificationGroup([cm], genome_build=self.genome_build) for cm in self.modifications]
        return None


class ClassificationGroups:

    def __init__(self,
                 classification_modifications: Iterable[ClassificationModification],
                 genome_build: Optional[GenomeBuild] = None):

        if not genome_build:
            genome_build = GenomeBuildManager.get_current_genome_build()
        self.genome_build = genome_build

        def clin_significance(cm: ClassificationModification) -> Optional[str]:
            return cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

        def condition_sorter(cm: ClassificationModification) -> Optional[str]:
            if resolved := cm.classification.condition_resolution_obj:
                if resolved.terms:
                    return "A" + (resolved.terms[0].name or resolved.terms[0].id).lower()
            return "Z" + (cm.get(SpecialEKeys.CONDITION) or "").lower()

        def condition_grouper(cm: ClassificationModification) -> Any:
            return cm.classification.condition_resolution_obj

        cached_chgvs: Dict[int, CHGVS] = dict()

        def c_hgvs_group(cm) -> CHGVS:
            if existing := cached_chgvs.get(cm.id):
                return existing
            c_hgvs = ClassificationGroup.c_hgvs_for(cm, genome_build)
            cached_chgvs[cm.id] = c_hgvs.without_transcript_version
            return c_hgvs

        evidence_keys: EvidenceKeyMap = EvidenceKeyMap.instance()

        groups: List[ClassificationGroup] = []

        # clinical significance, clin grouping, org
        sorted_by_clin_sig = list(classification_modifications)
        e_key_clin_sig = evidence_keys.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        sorted_by_clin_sig.sort(key=e_key_clin_sig.classification_sorter)
        for clin_sig, group1 in groupby(sorted_by_clin_sig, clin_significance):
            group1 = list(group1)
            # breakup by clinical grouping
            group1.sort(key=lambda cm: cm.classification.clinical_grouping_name)
            for _, group2 in groupby(group1, lambda cm: cm.classification.clinical_grouping_name):
                group2 = list(group2)
                # break up by org (TODO breakup by lab with optional breakup by org)
                group2.sort(key=lambda cm: cm.classification.lab.name)
                for _, group3 in groupby(group2, lambda cm: cm.classification.lab.name):
                    group3 = list(group3)
                    # breakup by transcript
                    group3.sort(key=c_hgvs_group)
                    for _, group4 in groupby(group3, c_hgvs_group):
                        group4 = list(group4)
                        group4.sort(key=condition_sorter)
                        for _, group5 in groupby(group4, condition_grouper):
                            actual_group = ClassificationGroup(modifications=group5, genome_build=genome_build, group_id=len(groups) + 1)
                            actual_group.clinical_significance_score = e_key_clin_sig.classification_sorter_value(clin_sig)
                            groups.append(actual_group)
        self.groups = groups
        self.groups.sort(key=lambda cg:cg.most_recent_curated)
        next_sort_order = 1
        for group in self.groups:
            group.sort_order = next_sort_order
            next_sort_order += 1
        self.groups.reverse()  # default store records as most recent to least recent

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
