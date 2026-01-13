from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from typing import TypedDict, Optional, Self

from django.db.models import TextChoices

from classification.criteria_strengths import CriteriaStrength
from classification.enums import AlleleOriginBucket, SpecialEKeys, CriteriaEvaluation, TestingContextBucket
from library.utils import strip_json

"""
return {
    "criteria_labels": self.criteria_labels,
    "classification_value": self.classification_value,
    "classification_sort": self.classification_sort,
    "pathogenicity": {
        "classification": self.classification_value,
        "sort": self.classification_sort,
    },
    "somatic": {
        "clinical_significance": self.somatic_clinical_significance,
        "amp_level": self.somatic_amp_level,
        "sort": self.somatic_sort
    },
    "date": {
        "value": curated_date.date_str,
        "type": curated_date.name
    }
}
"""


class ClassificationSummaryCacheDictPathogenicity(TypedDict):
    classification: Optional[str]
    sort: int
    # pending: Optional[str]


class ClassificationSummaryCacheDictSomatic(TypedDict):
    testing_context_bucket: Optional[str]
    tumor_type_category: Optional[str]  # condition grouping if testing_contest is solid-tumor
    clinical_significance: Optional[str]
    amp_level: Optional[str]
    sort: int


class ClassificationSummaryCachedDictDate(TypedDict):
    date: str
    type: Optional[str]


class ClassificationSummaryCacheDict(TypedDict):
    criteria_labels: list[str]
    pathogenicity: ClassificationSummaryCacheDictPathogenicity
    somatic: ClassificationSummaryCacheDictSomatic
    allele_origin_bucket: str
    date: ClassificationSummaryCachedDictDate


@dataclass(frozen=True)
class SomaticClinicalSignificanceValue:
    tier_level: str
    amp_level: Optional[str] = None

    @property
    def without_amp_level(self) -> Self:
        return SomaticClinicalSignificanceValue(tier_level=self.tier_level)

    @property
    def sort_value(self) -> Optional[int]:
        if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self):
            return sort_value
        elif self.amp_level:
            if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self.without_amp_level):
                return sort_value
        # default to 1, so things which shouldn't even be considered somatic can be 0
        return 1


_SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES = {
    SomaticClinicalSignificanceValue("tier_1", "A"): 10,
    SomaticClinicalSignificanceValue("tier_1", "B"): 9,
    SomaticClinicalSignificanceValue("tier_1"): 8,
    SomaticClinicalSignificanceValue("tier_1_or_2"): 7,
    SomaticClinicalSignificanceValue("tier_2", "C"): 6,
    SomaticClinicalSignificanceValue("tier_2", "D"): 5,
    SomaticClinicalSignificanceValue("tier_2"): 4,
    SomaticClinicalSignificanceValue("tier_3"): 3,
    SomaticClinicalSignificanceValue("tier_4"): 2
}


class ClassificationSummaryCalculator:

    def __init__(self, cm: 'ClassificationModification'):
        self.cm = cm

    def cache_dict(self) -> ClassificationSummaryCacheDict:
        from classification.models import CuratedDate
        curated_date = CuratedDate(self.cm).relevant_date

        pathology: ClassificationSummaryCacheDictPathogenicity = {
            "classification": self.classification_value,
            "sort": self.classification_sort,
            # "bucket": self.germline_bucket,
            # "pending": self.pending_classification_value
        }
        somatic: ClassificationSummaryCacheDictSomatic = {
            "testing_context_bucket": self.testing_context_bucket,
            "tumor_type_category": self.tumor_type_category,
            "clinical_significance": self.somatic_clinical_significance,
            "amp_level": self.somatic_amp_level,
            "sort": self.somatic_sort
        }
        date_json: ClassificationSummaryCachedDictDate = {
            "value": curated_date.date_str,
            "type": curated_date.name,
        }

        full_json: ClassificationSummaryCacheDict = {
            "criteria_labels": self.criteria_labels,
            "pathogenicity": pathology,
            "allele_origin_bucket": self.allele_origin_bucket,
            "somatic": somatic,
            "date": date_json
        }

        # strip out None values as that makes sorting work more naturally
        return strip_json(full_json)

    # @cached_property
    # def pending_classification_value(self) -> Optional[str]:
    #     from classification.models import classification_flag_types, ClassificationFlagTypes
    #     from flags.models import Flag, FlagStatus
    #     if flag := Flag.objects.filter(
    #         flag_type=classification_flag_types.classification_pending_changes,
    #         resolution__status=FlagStatus.OPEN,
    #         collection_id=self.cm.classification.flag_collection_id
    #     ).first():
    #         return flag.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY) if flag.data else 'Unknown'
    #     return None

    # @cached_property
    # def germline_bucket(self) -> Optional[int]:
    #     from classification.models import EvidenceKeyMap
    #     return EvidenceKeyMap.clinical_significance_to_bucket().get(self.classification_value)

    @cached_property
    def allele_origin_bucket(self) -> AlleleOriginBucket:
        return self.cm.classification.allele_origin_bucket

    @cached_property
    def is_possibly_somatic(self) -> bool:
        # TODO grab allele origin as it is per the modification
        return self.allele_origin_bucket != AlleleOriginBucket.GERMLINE

    @cached_property
    def somatic_amp_level(self) -> Optional[str]:
        if self.is_possibly_somatic and self.somatic_clinical_significance:
            for key, level in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
                if self.cm.get(key):
                    return level
        return None

    @cached_property
    def classification_value(self) -> Optional[str]:
        return self.cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @cached_property
    def classification_sort(self) -> int:
        if classification_value := self.classification_value:
            from classification.models import EvidenceKeyMap
            return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("sort_order").get(classification_value, 0)
        else:
            # sort unclassified to the end
            return None

    @cached_property
    def somatic_clinical_significance(self) -> Optional[str]:
        if self.is_possibly_somatic:
            return self.cm.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
        return None

    @cached_property
    def testing_context_bucket(self) -> TestingContextBucket:
        # TODO move this information into the testing context evidence key data
        if self.is_possibly_somatic:
            if testing_context_value := self.cm.get(SpecialEKeys.TESTING_CONTEXT):
                from classification.models import EvidenceKeyMap
                return EvidenceKeyMap.cached_key(SpecialEKeys.TESTING_CONTEXT).option_dictionary_property("testing_context_bucket").get(testing_context_value) or TestingContextBucket.OTHER
            else:
                return TestingContextBucket.UNKNOWN.value
        else:
            return TestingContextBucket.GERMLINE.value

    @cached_property
    def tumor_type_category(self) -> Optional[str]:
        if self.testing_context_bucket == TestingContextBucket.SOLID_TUMOR:
            # TODO actually calculate
            return None
        elif self.testing_context_bucket == TestingContextBucket.HAEMATOLOGY:
            # TODO actually calculate
            return None
        return None

    @cached_property
    def somatic_sort(self) -> Optional[int]:
        if self.is_possibly_somatic:
            return SomaticClinicalSignificanceValue(
                self.somatic_clinical_significance,
                self.somatic_amp_level
            ).sort_value
        else:
            return None

    @cached_property
    def criteria_labels(self) -> list[str]:
        from classification.models import EvidenceKeyMap
        cm = self.cm
        strengths: set[CriteriaStrength] = set()
        for e_key in EvidenceKeyMap.cached().criteria():
            strength = cm.get(e_key.key)
            if CriteriaEvaluation.is_met(strength):
                strengths.add(CriteriaStrength(e_key, strength))
        for amp_level, letter in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
            if value := cm.get_value_list(amp_level):
                e_key = EvidenceKeyMap.cached_key(amp_level)
                for sub_value in value:
                    sub_value_label = e_key.pretty_value(sub_value)
                    strengths.add(CriteriaStrength(
                        ekey=EvidenceKeyMap.cached_key(amp_level),
                        custom_strength=f"{letter}_{sub_value_label}")
                    )
        return list(str(x) for x in sorted(strengths))
