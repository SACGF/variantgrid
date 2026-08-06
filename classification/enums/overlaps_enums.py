from dataclasses import dataclass
from typing import Optional
from dataclasses_json import DataClassJsonMixin
from django.db.models import TextChoices
from django.utils.safestring import mark_safe

from classification.enums import SpecialEKeys, OverlapStatus, OverlapOverrideStatus


class OverlapType(TextChoices):
    SINGLE_CONTEXT = "context", "Single Context"
    CROSS_CONTEXT = "cross", "Cross Context"

    @property
    def priority_order(self) -> int:
        match self:
            case OverlapType.SINGLE_CONTEXT: return 1
            case OverlapType.CROSS_CONTEXT: return 2

    def __lt__(self, other):
        return self.priority_order < other.priority_order


class ClassificationResultValue(TextChoices):
    """
    Represents the different outcome values types of a Classification
    e.g. a Germline or Somatic record will have a Pathogenicity or Oncogenicity score (comparable so treated as the same value type)
    and a Somatic record can also have
    """

    ONC_PATH = "O", "Onco-Path"
    SOMATIC_CLINICAL_SIGNIFICANCE = "S", "Clinical significance"

    @staticmethod
    @property
    def supported_properties():
        from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
        if not OVERLAP_CLIN_SIG_ENABLED:
            return [ClassificationResultValue.ONC_PATH]
        else:
            return [ClassificationResultValue.ONC_PATH, ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE]

    @property
    def priority_order(self) -> int:
        match self:
            case ClassificationResultValue.ONC_PATH: return 1
            case ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE: return 2

    @property
    def evidence_key_str(self) -> str:
        match self:
            case ClassificationResultValue.ONC_PATH: return SpecialEKeys.ONC_PATH
            case ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE: return SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE


class OverlapContributionStatus(TextChoices):
    """
    Represents if an OverlapContribution (typically a ClassificationGrouping) is contributing to the Overlap.
    """

    PENDING_CALCULATION = "P", "Pending Calculation"
    CONTRIBUTING = "C", "Contributing"
    NOT_SHARED = "N", "Not-shared"
    NO_VALUE = "X", "No value"
    NON_COMPARABLE_VALUE = "Z", "Non-comparable value"  # e.g. Risk Factor


class TriageStatus(TextChoices):
    """
    An OverlapContribution has a triage status, generally saying how confident they are in the record
    """

    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    AMENDED = "A", "Amended"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Value"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"
    NON_INTERACTIVE_THIRD_PARTY = "Z", "Non-Interactive Party"

    @property
    def is_default(self):
        # default as in the user likely hasn't changed anything
        return self in {TriageStatus.PENDING, TriageStatus.NON_INTERACTIVE_THIRD_PARTY}

    @property
    def icon(self):
        match self:
            case TriageStatus.PENDING:
                return mark_safe('<i class="fa-regular fa-clock"  title="Pending" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_WILL_FIX:
                return mark_safe('<i class="fa-solid fa-strikethrough" title="Will Amend" style="opacity:0.6"></i>')
            case TriageStatus.AMENDED:
                return mark_safe('<i class="fa-solid fa-strikethrough" title="Have Amended" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_WILL_DISCUSS:
                return mark_safe('<i class="fa-regular fa-comments" title="For Joint Discussion" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_SATISFACTORY:
                return mark_safe('<i class="fa-solid fa-clipboard-check" title="Confident in Value" style="opacity:0.6"></i>')
            case TriageStatus.COMPLEX:
                return mark_safe('<i class="fa-solid fa-clipboard-question" title="Complex Reasons for Discordance" style="opacity:0.6"></i>')
            case TriageStatus.NON_INTERACTIVE_THIRD_PARTY:
                return mark_safe('<i class="fa-solid fa-face-meh-blank" title="3rd Party" style="opacity:0.6"></i>')
            case _:
                return "?"


class OverlapEntrySourceTextChoices(TextChoices):
    CLASSIFICATION = "CLASS", "CLASSIFICATION"
    CLINVAR = "CLIN", "CLINVAR"


@dataclass(frozen=True)
class TriageState(DataClassJsonMixin):
    """
    Combination of a TriageStatus, and an optional amend_value (only if TriageStatus is REVIEWED_WILL_FIX)
    e.g. PENDING or REVIEWED WILL FIX - (P)athogenic
    """
    status: TriageStatus = TriageStatus.PENDING
    amend_value: Optional[str] = None

    def __str__(self):
        if self.amend_value:
            # TODO does amend_value need to be formatted?
            return f"{self.status.label} ({self.amend_value.replace("_", "-")})"
        return self.status.label

    @staticmethod
    def default_json():
        return TriageState().to_dict()


@dataclass(frozen=True)
class TriageComment(DataClassJsonMixin):
    """
    A comment plus a counter of how many comments there's been
    """
    text: Optional[str] = None
    count: int = 0

    def next_comment(self, text: Optional[str] = None):
        return TriageComment(
            text=text,
            count=self.count + 1,
        )

    @staticmethod
    def default_json():
        return TriageComment().to_dict()

    def __bool__(self):
        return bool(self.text)


@dataclass(frozen=True)
class OverlapState(DataClassJsonMixin):
    status: OverlapStatus
    has_pending_values: bool = False
    override_status: OverlapOverrideStatus = OverlapOverrideStatus.NO_OVERRIDE

    @property
    def label(self):
        if override_status := self.override_status:
            return f"Reviewed as {override_status.label}"
        else:
            return self.status.label

    @property
    def is_active_discordance(self):
        return self.status.is_discordant and not self.override_status

    @staticmethod
    def default_json():
        return OverlapState(status=OverlapStatus.NO_CONTRIBUTIONS, has_pending_values=False, override_status=OverlapOverrideStatus.NO_OVERRIDE).to_dict()
