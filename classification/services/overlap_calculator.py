from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Iterable
from classification.enums import OverlapStatus, OverlapOverrideStatus, TriageStatus
from classification.models import ClassificationResultValue, \
    EvidenceKeyMap, OverlapContribution, ClassificationSummaryCacheObj
from classification.enums.overlaps_enums import OverlapContributionStatus
from library.utils import first


OVERLAP_CLIN_SIG_ENABLED = False  # ensure reference this whenever doing some functionality when ClinSign should be supported
# So it's easy to work out where to write functionality when we start supporting it
# The idea is NOT to have some environments support it and some not


@dataclass(frozen=True)
class OverlapStatusCalculation:
    current_value: OverlapStatus
    has_pending_values: bool = False
    override_status: Optional[OverlapOverrideStatus] = OverlapOverrideStatus.NO_OVERRIDE

    @property
    def effective_value(self):
        return self.override_pending_value if (self.override_pending_value is not None) else self.current_value


class OverlapCalculatorBase(ABC):

    @classmethod
    @abstractmethod
    def value_type(cls) -> ClassificationResultValue:
        raise NotImplementedError()

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheObj) -> Optional[str]:
        raise NotImplementedError()

    @classmethod
    def is_comparable_value(cls, value):
        return True

    @classmethod
    def calculate_entries(cls, entries: Iterable[OverlapContribution]) -> OverlapStatusCalculation:
        non_comparable_values: int = 0
        contributing: list[OverlapContribution] = []

        for entry in entries:
            # FIXME, looking for a different solution that possibly outdated
            # TODO do we want to do a cleaner solution than this for ExpertPanels?
            # if entry.possibly_outdated:
            #     continue

            match entry.contribution_status:
                case OverlapContributionStatus.CONTRIBUTING:
                    contributing.append(entry)
                case OverlapContributionStatus.NON_COMPARABLE_VALUE:
                    non_comparable_values += 1
                case _:
                    pass  # don't care about unshared for calculation values

        has_pending_values = any(con.is_amending for con in contributing)
        if len(contributing) == 0:
            if non_comparable_values > 0:
                return OverlapStatusCalculation(OverlapStatus.NO_COUNTING_CONTRIBUTIONS, has_pending_values=has_pending_values)
            else:
                return OverlapStatusCalculation(OverlapStatus.NO_CONTRIBUTIONS, has_pending_values=has_pending_values)
        elif len(contributing) == 1:
            return OverlapStatusCalculation(OverlapStatus.SINGLE_SUBMITTER, has_pending_values=has_pending_values)
        else:
            override_value: OverlapOverrideStatus = OverlapOverrideStatus.NO_OVERRIDE

            all_values = set(con.effective_value for con in contributing)
            base_value: OverlapStatus
            if len(all_values) == 1:
                base_value = OverlapStatus.EXACT_AGREEMENT
            else:
                base_value = cls.calculate_status_for_multiple_entries(all_values)

            third_party = [con for con in contributing if con.triage_state_obj.status == TriageStatus.NON_INTERACTIVE_THIRD_PARTY]
            interactive_contributors = [con for con in contributing if con.triage_state_obj.status != TriageStatus.NON_INTERACTIVE_THIRD_PARTY]

            if interactive_contributors:
                all_reviewed = all(con.is_review_agreed_value_met() for con in interactive_contributors)
                all_complex = all(con.triage_state_obj.status == TriageStatus.COMPLEX for con in interactive_contributors)

                if all_reviewed:
                    override_value = OverlapOverrideStatus.CONTINUED_DISCORDANCE
                elif all_complex:
                    override_value = OverlapOverrideStatus.COMPLEX
                elif base_value.is_discordant:
                    # see if it's ClinVar that's making the over discordant
                    non_clinvar_values = set(con.effective_value for con in interactive_contributors)
                    clinvar_discordant = len(non_clinvar_values) > 1 and cls.calculate_status_for_multiple_entries(all_values).is_discordant

                    if clinvar_discordant:
                        max_clinvar_date = max(con.effecive_date for con in third_party)
                        max_classification_date = max(con.effecive_date for con in interactive_contributors)
                        if max_clinvar_date < max_classification_date:
                            override_value = OverlapOverrideStatus.IGNORING_OLD_CLINVAR
                        elif all(con.triage_state_obj.status == TriageStatus.REVIEWED_SATISFACTORY for con in interactive_contributors): # all confident
                            override_value = OverlapOverrideStatus.CONFIDENT_VS_CLINVAR

            return OverlapStatusCalculation(base_value, has_pending_values, override_value)

    @classmethod
    @abstractmethod
    def calculate_status_for_multiple_entries(cls, values: set[str]) -> OverlapStatus:
        raise NotImplementedError()


class OverlapCalculatorClinSig(OverlapCalculatorBase):

    @classmethod
    def value_type(cls) -> ClassificationResultValue:
        return ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheObj) -> Optional[str]:
        return summary.somatic.clinical_significance

    @classmethod
    def is_comparable_value(cls, value):
        return True

    @classmethod
    def calculate_status_for_multiple_entries(cls, values: set[str]) -> OverlapStatus:
        has_tier_1_and_2 = False
        tiers = set()
        for value in values:
            if value == "tier_1_or_2":
                has_tier_1_and_2 = True
            else:
                tiers.add(value)

        if tiers == {"tier_1", "tier_2"}:
            return OverlapStatus.TIER_1_VS_TIER_2_DIFFERENCES
        elif len(tiers) == 1 and first(tiers) in {"tier_1", "tier_2"} and has_tier_1_and_2:
            return OverlapStatus.RESOLUTION_DIFFERENCES
        else:
            return OverlapStatus.MEDICALLY_SIGNIFICANT


class OverlapCalculatorOncPath(OverlapCalculatorBase):

    @classmethod
    def value_type(cls) -> ClassificationResultValue:
        return ClassificationResultValue.ONC_PATH

    @classmethod
    def is_comparable_value(cls, value):
        if EvidenceKeyMap.clinical_significance_to_bucket().get(value) is None:
            return False
        return True

    @classmethod
    def value_from_summary(cls, summary: ClassificationSummaryCacheObj) -> Optional[str]:
        return summary.pathogenicity.classification

    @classmethod
    def calculate_status_for_multiple_entries(cls, values: set[str]) -> OverlapStatus:
        """
        :param values: 2+ OverlapEntries all contributing, should have at least 1 difference
        :return: The calculated Overlap Status for Onc or Pathogenicity
        """

        # must be 2 or more entries, all entries should be CONTRIBUTING
        all_classification_values: set[str] = set()
        all_bucket_values: set[int] = set()

        for value in values:
            bucket = EvidenceKeyMap.clinical_significance_to_bucket().get(value)
            all_bucket_values.add(bucket)
            all_classification_values.add(value)

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
    elif value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
        return OverlapCalculatorClinSig()
    else:
        raise ValueError(f"Unsupported value type {value_type}")
