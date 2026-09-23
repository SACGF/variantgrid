from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Self
from dataclasses_json import DataClassJsonMixin
from classification.criteria_strengths import CriteriaStrength
from classification.enums import AlleleOriginBucket, SpecialEKeys, CriteriaEvaluation, TestingContextBucket, ClassificationResultValue


@dataclass(frozen=True)
class ClassificationSummaryCacheObjDate(DataClassJsonMixin):
    # TODO, deprecate in favour of EffectiveDate
    date: str  # in format of yyyy-mm-dd
    type: Optional[str]

    def __bool__(self):
        return bool(self.date)

    def __lt__(self, other):
        return self.date < other.date


@dataclass(frozen=True)
class ClassificationSummaryCacheObjSomatic(DataClassJsonMixin):
    testing_context_bucket: Optional[str]
    tumor_type_category: Optional[str]  # condition grouping if testing_contest is solid-tumor
    clinical_significance: Optional[str]
    amp_level: Optional[str]
    sort: Optional[int]

    @property
    def somatic_clinical_significance_value(self) -> Optional['SomaticClinicalSignificanceValue']:
        if cs := self.clinical_significance:
            return SomaticClinicalSignificanceValue(cs, self.amp_level)
        return None


@dataclass(frozen=True)
class ClassificationSummaryCacheObjPathogenicity(DataClassJsonMixin):
    classification: Optional[str]
    sort: Optional[int]

    @property
    def bucket(self):
        from classification.models import EvidenceKeyMap
        return EvidenceKeyMap.onc_path_to_bucket().get(self.classification)


@dataclass(frozen=True)
class ClassificationSummaryCacheObj(DataClassJsonMixin):
    criteria_labels: list[str]
    pathogenicity: ClassificationSummaryCacheObjPathogenicity
    somatic: ClassificationSummaryCacheObjSomatic
    allele_origin_bucket: str
    date: ClassificationSummaryCacheObjDate

    def value_for_value_type(self, value_type: ClassificationResultValue):
        match value_type:
            case ClassificationResultValue.ONC_PATH: return self.pathogenicity.classification
            case ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE: return self.somatic.clinical_significance
            case _: raise ValueError(f"Unexpected value type {value_type}")

    @staticmethod
    def from_dict_safe(data: dict):
        # is there a way we can do this with schema?
        # the problem is we don't want to substitute with None, but empty objects
        data_copy = dict(data)
        for man_sub_field in ["pathogenicity", "somatic", "date"]:
            if man_sub_field not in data_copy:
                data_copy[man_sub_field] = {}
        return ClassificationSummaryCacheObj.from_dict(data_copy, infer_missing=True)


@dataclass(frozen=True)
class SomaticClinicalSignificanceValue:
    tier_level: str
    amp_level: Optional[str] = None

    @property
    def without_amp_level(self) -> Self:
        return SomaticClinicalSignificanceValue(tier_level=self.tier_level)

    @property
    def sort_value(self) -> int:
        if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self):
            return sort_value
        elif self.amp_level:
            if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self.without_amp_level):
                return sort_value
        # default to 1, so things which shouldn't even be considered somatic can be 0
        return 1

    @property
    def pretty_value(self):
        from classification.models import EvidenceKeyMap
        parts = [EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(self.tier_level)]
        if amp_level := self.amp_level:
            parts.append(amp_level)
        return "".join(parts)

    def __lt__(self, other):
        return self.sort_value < other.sort_value


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

    def cache_dict(self) -> dict:
        from classification.models import CuratedDate
        curated_date = CuratedDate(self.cm).relevant_date

        full_obj = ClassificationSummaryCacheObj(
            criteria_labels=self.criteria_labels,
            pathogenicity=ClassificationSummaryCacheObjPathogenicity(
                classification=self.classification_value,
                sort=self.classification_sort
            ),
            somatic=ClassificationSummaryCacheObjSomatic(
                testing_context_bucket=self.testing_context_bucket,
                tumor_type_category=self.tumor_type_category,
                clinical_significance=self.somatic_clinical_significance,
                amp_level=self.somatic_amp_level,
                sort=self.somatic_sort
            ),
            allele_origin_bucket=self.allele_origin_bucket,
            date=ClassificationSummaryCacheObjDate(
                date=curated_date.date_str,
                type=curated_date.name
            )
        )
        return full_obj.to_dict()

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


# FIXME underlying code should change due to overlaps
def clinical_significance_pills(summary: ClassificationSummaryCacheObj, allele_origin_bucket: str) -> list[dict]:
    from classification.models import EvidenceKeyMap
    """ Label/CSS class for the c-pill spans - rendered server side by the clinical_significance_values
        tag, and client side by the variant details samples grid. CSS is .c-pill.cs-* / .c-pill.scs-* """
    pathogenicity = summary.pathogenicity
    somatic = summary.somatic

    germline_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    pending_from = None
    value = pathogenicity.classification
    # TODO pending from isn't in the ClassificationSummary anymore

    pills = [{
        "title": germline_key.pretty_label,
        "pending_from": pending_from,
        "label": germline_key.pretty_value(value, value) or "No Data",
        "css_class": "cs cs-" + (value.lower() if value else "none"),
    }]

    always_show_somatic = allele_origin_bucket != AlleleOriginBucket.GERMLINE
    if always_show_somatic or somatic.clinical_significance:
        somatic_key = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
        value = somatic.clinical_significance
        label = somatic_key.pretty_value(value, value) or "No Data"
        if amp_level := somatic.amp_level:
            label += amp_level
        pills.append({
            "title": somatic_key.pretty_label,
            "label": label,
            "css_class": f"scs scs-{value.lower() if value else 'none'}",
        })
    return pills
