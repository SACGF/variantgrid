"""
ClinicalContext vocabulary: ClinicalContextStatus, DiscordanceLevel (how far apart the shared
classifications in a context are) and ClinicalContextRecalcTrigger (what caused a recalculation).
"""
from enum import Enum


class ClinicalContextStatus:
    CONCORDANT = 'C'
    DISCORDANT = 'D'

    CHOICES = (
        (CONCORDANT, 'Concordant'),
        (DISCORDANT, '*Discordant*')
    )


class ClinicalContextRecalcTrigger(Enum):
    ADMIN = "admin"
    VARIANT_SET = "variant_set"
    CLINICAL_GROUPING_SET = "clinical_grouping"
    WITHDRAW_DELETE = "withdraw_delete"
    SUBMISSION = "submission"
    DELAYED = "delayed"
    OTHER = "other"
    PENDING_CS_CHANGE = "pending_cs_change"
    RE_OPEN = "reopen"
    CONTINUED_DISCORDANCE = "continued_discordance"


class DiscordanceLevel(str, Enum):
    """
    Values are assigned based on how big the differences are considered when it comes to discordance
    e.g. difference between Benign and Likely Benign isn't as important as the difference between Likely Benign and VUS etc
    """

    NO_ENTRIES = 'no_entries'
    SINGLE_SUBMISSION = 'single_submission'  # single shared submission
    CONCORDANT_AGREEMENT = 'concordant_agreement'  # complete agreement
    CONCORDANT_DIFF_VUS = 'concordant_agreement_diff_vus'  # VUS-A, VUS-B vs VUS-C
    CONCORDANT_CONFIDENCE = 'concordant_confidence'  # Benign vs Likely Benign
    MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED = 'multiple_records'  # Somatic before we support discordan
    DISCORDANT = 'discordant'

    @property
    def label(self) -> str:
        if self == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return "Concordant (Agreement)"
        if self == DiscordanceLevel.CONCORDANT_DIFF_VUS:
            return "Concordant (Agreement Differing VUS)"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "Concordant (Confidence)"
        if self == DiscordanceLevel.NO_ENTRIES:
            return "No Shared Submissions"
        if self == DiscordanceLevel.SINGLE_SUBMISSION:
            return "Single Shared Submission"
        if self == DiscordanceLevel.DISCORDANT:
            return "Discordant"
        if self == DiscordanceLevel.MULTIPLE_RECORDS_DISCORDANCE_NOT_SUPPORTED:
            return "Multiple Records"
        return "Unknown"

    @property
    def css_class(self) -> str:
        return f"overlap-{self.value}"
