from django.db.models import TextChoices
from django.utils.safestring import mark_safe


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
    REVIEWED_SATISFACTORY = "R", "Confident in Classification"
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
            case _:
                return "?"


class OverlapEntrySourceTextChoices(TextChoices):
    CLASSIFICATION = "CLASS", "CLASSIFICATION"
    CLINVAR = "CLIN", "CLINVAR"
