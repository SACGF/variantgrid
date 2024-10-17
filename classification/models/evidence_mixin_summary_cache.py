from dataclasses import dataclass
from functools import cached_property
from typing import TypedDict, Optional, Self
from classification.criteria_strengths import CriteriaStrength
from classification.enums import AlleleOriginBucket, SpecialEKeys, CriteriaEvaluation

"""
return {
            "criteria_labels": self.criteria_labels,
            "classification_value": self.classification_value,
            "classification_sort": self.classification_sort,
            "pathogenicity": {
                "classification": self.classification_value,
                "sort": self.classification_sort,
                "bucket": self.germline_bucket
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
    bucket: Optional[int]


class ClassificationSummaryCacheDictSomatic(TypedDict):
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
        return None


_SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES = {
    SomaticClinicalSignificanceValue("tier_1", "A"): 9,
    SomaticClinicalSignificanceValue("tier_1", "B"): 8,
    SomaticClinicalSignificanceValue("tier_1"): 7,
    SomaticClinicalSignificanceValue("tier_1_or_2"): 6,
    SomaticClinicalSignificanceValue("tier_2", "C"): 5,
    SomaticClinicalSignificanceValue("tier_2", "D"): 4,
    SomaticClinicalSignificanceValue("tier_2"): 3,
    SomaticClinicalSignificanceValue("tier_3"): 2,
    SomaticClinicalSignificanceValue("tier_4"): 1
}


class ClassificationSummaryCalculator:

    def __init__(self, cm: 'ClassificationModification'):
        self.cm = cm

    def cache_dict(self) -> ClassificationSummaryCacheDict:
        from classification.models import CuratedDate
        curated_date = CuratedDate(self.cm).relevant_date

        return {
            "criteria_labels": self.criteria_labels,
            "classification_value": self.classification_value,
            "classification_sort": self.classification_sort,
            "pathogenicity": {
                "classification": self.classification_value,
                "sort": self.classification_sort,
                "bucket": self.germline_bucket
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

    @cached_property
    def germline_bucket(self) -> Optional[int]:
        from classification.models import EvidenceKeyMap
        return EvidenceKeyMap.clinical_significance_to_bucket().get(self.classification_value)

    @cached_property
    def is_possibly_somatic(self) -> bool:
        return self.cm.classification.allele_origin_bucket != AlleleOriginBucket.GERMLINE

    @cached_property
    def somatic_amp_level(self) -> Optional[str]:
        if self.is_possibly_somatic:
            for key, level in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
                if self.cm.get(key):
                    return level

    @cached_property
    def classification_value(self) -> Optional[str]:
        return self.cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @cached_property
    def classification_sort(self) -> int:
        if classification_value := self.classification_value:
            from classification.models import EvidenceKeyMap
            return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_indexes.get(classification_value, 0) + 1
        else:
            return 0

    @cached_property
    def somatic_clinical_significance(self) -> Optional[str]:
        if self.is_possibly_somatic:
            return self.cm.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)

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