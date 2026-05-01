from dataclasses import dataclass
from datetime import datetime, date
from typing import Optional, Union

from dataclasses_json import DataClassJsonMixin
from django.db.models import TextChoices
from django.utils.safestring import mark_safe

from classification.models import CuratedDate



class OverlapType(TextChoices):
    SINGLE_CONTEXT = "context", "Single Context"
    CROSS_CONTEXT = "cross", "Cross Context"
    # CLINVAR_EXPERT_PANEL = "clinvar", "ClinVar Expert Panel"
    # CROSS_GERMLINE_NON_CANCER = "cross_germ_non_cancer", "Cross Germline / Somatic Non-Cancer"
    # CROSS_GERMLINE_HAEM = "cross_germ_haem", "Cross Germline / Somatic Haem"
    # CROSS_GERMLINE_SOLID_TUMOR = "cross_germ_solid_tumor", "Cross Germline / Solid Tumour"
    # CROSS_HAEM_SOLID_TUMOR = "cross_haem_solid_tumor", "Cross Haem / Solid Tumour"

    @property
    def priority_order(self) -> int:
        match self:
            case OverlapType.SINGLE_CONTEXT: return 1
            # case OverlapType.CLINVAR_EXPERT_PANEL: return 2
            case OverlapType.CROSS_CONTEXT: return 3

    def __lt__(self, other):
        return self.priority_order < other.priority_order


class ClassificationResultValue(TextChoices):
    # FIXME should be called value type
    ONC_PATH = "O", "Onco-Path"
    CLINICAL_SIGNIFICANCE = "S", "Clinical significance"

    @staticmethod
    @property
    def supported_properties():
        from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
        if not OVERLAP_CLIN_SIG_ENABLED:
            return [ClassificationResultValue.ONC_PATH]
        else:
            return [ClassificationResultValue.ONC_PATH, ClassificationResultValue.CLINICAL_SIGNIFICANCE]

    @property
    def priority_order(self) -> int:
        match self:
            case ClassificationResultValue.ONC_PATH: return 1
            case ClassificationResultValue.CLINICAL_SIGNIFICANCE: return 2


class OverlapContributionStatus(TextChoices):
    PENDING_CALCULATION = "P", "Pending Calculation"
    CONTRIBUTING = "C", "Contributing"
    NOT_SHARED = "N", "Not-shared"
    NO_VALUE = "X", "No value"
    NON_COMPARABLE_VALUE = "Z", "Non-comparable value"  # e.g. Risk Factor


class TriageStatus(TextChoices):
    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Value"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"
    NON_INTERACTIVE_THIRD_PARTY = "Z", "Non-Interactive Party"

    @property
    def icon(self):
        match self:
            case TriageStatus.PENDING:
                return mark_safe('<i class="fa-regular fa-clock"  title="Pending" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_WILL_FIX:
                return mark_safe('<i class="fa-solid fa-strikethrough" title="Will Amend" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_WILL_DISCUSS:
                return mark_safe('<i class="fa-regular fa-comments" title="For Joint Discussion" style="opacity:0.6"></i>')
            case TriageStatus.REVIEWED_SATISFACTORY:
                return mark_safe('<i class="fa-solid fa-clipboard-check" title="Confident in Classification" style="opacity:0.6"></i>')
            case TriageStatus.COMPLEX:
                return mark_safe('<i class="fa-solid fa-clipboard-question" title="Complex Reasons for Discordance" style="opacity:0.6"></i>')
            case TriageStatus.NON_INTERACTIVE_THIRD_PARTY:
                return mark_safe('<i class="fa-solid fa-face-meh-blank" title="3rd Party" style="opacity:0.6"></i>')
            case _:
                return "?"


class OverlapEntrySourceTextChoices(TextChoices):
    CLASSIFICATION = "CLASS", "CLASSIFICATION"
    CLINVAR = "CLIN", "CLINVAR"


@dataclass
class TriageState(DataClassJsonMixin):
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

    def __str__(self):
        return self.text