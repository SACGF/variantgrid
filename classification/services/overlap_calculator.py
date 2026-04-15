from abc import ABC, abstractmethod
from typing import Optional, Iterable
from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.models import ClinVarRecord
from annotation.templatetags.clinvar_tags import ClinVarDetails
from classification.enums import TestingContextBucket, OverlapStatus
from classification.models import ClassificationResultValue, ClassificationSummaryCacheDict, \
    EvidenceKeyMap, OverlapContribution, TriageStatus
from classification.models.overlaps_enums import OverlapContributionStatus, OverlapEntrySourceTextChoices
from library.utils import first
from snpdb.models import Allele


OVERLAP_CLIN_SIG_ENABLED = False  # ensure reference this whenever doing some functionality when ClinSign should be supported
# So it's easy to work out where to write functionality when we start supporting it
# The idea is NOT to have some environments support it and some not


class OverlapCalculatorBase(ABC):

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        raise NotImplementedError()

    @classmethod
    def clinvar_to_entry(cls, allele: Allele) -> Optional[OverlapContribution]:
        raise NotImplementedError()

    @classmethod
    def is_comparable_value(cls, value):
        return True

    @classmethod
    def calculate_entries(cls, entries: Iterable[OverlapContribution]) -> OverlapStatus:
        non_comparable_values: int = 0
        contributing: list[OverlapContribution] = []
        for entry in entries:
            match entry.contribution_status:
                case OverlapContributionStatus.CONTRIBUTING:
                    contributing.append(entry)
                case OverlapContributionStatus.NON_COMPARABLE_VALUE:
                    non_comparable_values += 1
                case _:
                    pass  # don't care about unshared for calculation values

        if len(contributing) == 0:
            if non_comparable_values > 0:
                return OverlapStatus.NO_COUNTING_CONTRIBUTIONS
            else:
                return OverlapStatus.NO_CONTRIBUTIONS
        elif len(contributing) == 1:
            return OverlapStatus.SINGLE_SUBMITTER
        else:
            if all(con.value == contributing[0].value for con in contributing):
                return OverlapStatus.EXACT_AGREEMENT
            else:
                return cls._calculate_status_for_multiple_entries(contributing)

    @classmethod
    @abstractmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        raise NotImplementedError()


class OverlapCalculatorClinSig(OverlapCalculatorBase):

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        return summary.get("somatic", {}).get("clinical_significance")

    @classmethod
    def is_comparable_value(cls, value):
        return True

    @classmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        has_tier_1_and_2 = False
        tiers = set()
        for entry in entries:
            if entry.value == "tier_1_or_2":
                has_tier_1_and_2 = True
            else:
                tiers.add(entry.value)

        if tiers == {"tier_1", "tier_2"}:
            return OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES
        elif len(tiers) == 1 and first(tiers) in {"tier_1", "tier_2"} and has_tier_1_and_2:
            return OverlapStatus.RESOLUTION_DIFFERENCES
        else:
            return OverlapStatus.MEDICALLY_SIGNIFICANT


class OverlapCalculatorOncPath(OverlapCalculatorBase):

    @classmethod
    def is_comparable_value(cls, value):
        if EvidenceKeyMap.clinical_significance_to_bucket().get(value) is None:
            return False
        return True

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheDict) -> str:
        return summary.get("pathogenicity", {}).get("classification")

    @classmethod
    def clinvar_to_contribution(cls, allele: Allele) -> Optional[OverlapContribution]:
        # FIXME not used
        if clinvar_details := ClinVarDetails.instance_from(allele=allele):
            if clinvar_details.is_expert_panel_or_greater and clinvar_details.clinvar.highest_pathogenicity > 0:
                clinvar_record_collection = ClinVarFetchRequest(
                    clinvar_variation_id=clinvar_details.clinvar.clinvar_variation_id,
                ).fetch()
                expert_panel: ClinVarRecord
                if expert_panel := clinvar_record_collection.expert_panel:
                    value = expert_panel.clinical_significance
                    relevant_value = ClassificationResultValue.ONC_PATH and EvidenceKeyMap.clinical_significance_to_bucket().get(value) is not None
                    effective_date = expert_panel.date_last_evaluated or expert_panel.date_clinvar_updated
                    print(effective_date)

                    oc = OverlapContribution(
                        source=OverlapEntrySourceTextChoices.CLINVAR,
                        scv=expert_panel.record_id,
                        testing_context_bucket=TestingContextBucket.GERMLINE,
                        allele=allele,
                        classification_grouping_id=None,
                        value=value,
                        contribution_status=OverlapContributionStatus.CONTRIBUTING if relevant_value else OverlapContributionStatus.NON_COMPARABLE_VALUE,
                        effective_date=effective_date,
                        # effective_date_type= # FIXME
                        triage_status=TriageStatus.NON_INTERACTIVE_THIRD_PARTY
                    )
                    return oc
        return None

    @classmethod
    def _calculate_status_for_multiple_entries(cls, entries: list[OverlapContribution]) -> OverlapStatus:
        """
        :param entries: 2+ OverlapEntries all contributing, should have at least 1 difference
        :return: The calculated Overlap Status for Onc or Pathogenicity
        """

        # must be 2 or more entries, all entries should be CONTRIBUTING
        all_classification_values: set[str] = set()
        all_bucket_values: set[int] = set()

        for entry in entries:
            bucket = EvidenceKeyMap.clinical_significance_to_bucket().get(entry.value)
            all_bucket_values.add(bucket)
            all_classification_values.add(entry.value)

        desired_value: OverlapStatus
        if len(all_classification_values) == 2 and "VUS" in all_classification_values and len(all_bucket_values) == 1:
            # here we would have VUS and one of VUS_A, VUS_B, VUS_C
            return OverlapStatus.RESOLUTION_DIFFERENCES
        elif all_classification_values == {"P", "O"} or all_classification_values == {"LP", "LO"}:
            return OverlapStatus.TERMINOLOGY_DIFFERENCES
        elif len(all_bucket_values) == 1:
            return OverlapStatus.MINOR_DIFFERENCES
        elif len(all_bucket_values) > 1:
            if 3 in all_bucket_values:
                return OverlapStatus.MEDICALLY_SIGNIFICANT
            else:
                return OverlapStatus.MAJOR_DIFFERENCES
        else:
            raise ValueError("Unhandled calculation state")


def calculator_for_value_type(value_type: ClassificationResultValue) -> OverlapCalculatorBase:
    if value_type == ClassificationResultValue.ONC_PATH:
        return OverlapCalculatorOncPath()
    elif value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
        return OverlapCalculatorClinSig()
    else:
        raise ValueError(f"Unsupported value type {value_type}")
