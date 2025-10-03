"""
Copying the Discordance Model Utils but for conflicts
"""
from dataclasses import dataclass
from functools import cached_property

from classification.enums import ConflictSeverity
from classification.models import Conflict, DiscordanceReportNextStep
from library.utils import first
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab


@dataclass
class ConflictCategoriesUtils:
    active: int = 0
    medical: int = 0
    waiting_for_triage: int = 0
    waiting_for_triage_medical: int = 0
    waiting_for_amend: int = 0
    ready_to_discuss: int = 0


class ConflictReportCategories:

    def __init__(self, perspective: LabPickerData):
        self.perspective = perspective

        # discordant_c = DiscordanceReportClassification.objects \
        #     .filter(classification_original__classification__lab__in=perspective.selected_labs) \
        #     .values_list('report_id', flat=True)
        # # .filter(classification_original__classification__withdrawn=False)  used to
        # self.dr_qs = DiscordanceReport.objects.filter(pk__in=discordant_c)\
        #     .order_by('-created')\
        #     .prefetch_related('discordancereporttriage_set')\
        #     .prefetch_related('discordancereportclassification_set')\
        #     .select_related('clinical_context')

    def labs(self) -> list[Lab]:
        return sorted(self.perspective.your_labs)

    def labs_quick_str(self) -> str:
        if len(self.perspective.selected_labs) == 1:
            return str(first(self.perspective.selected_labs))
        else:
            return "your assigned labs"

    @cached_property
    def all_counts(self) -> ConflictCategoriesUtils:
        lab_ids = LabPickerData.for_user(self.user, self.get_query_param("lab")).lab_ids
        qs = Conflict.objects.all().for_labs(self.perspective.selected_labs)

        # TODO need to distinguish medically significant
        qs = qs.filter(severity__gte=ConflictSeverity.MAJOR)

        # UNANIMOUSLY_COMPLEX = "C"
        # AWAITING_YOUR_TRIAGE = "T"
        # AWAITING_YOUR_AMEND = "A"
        # AWAITING_OTHER_LAB = "O"
        # TO_DISCUSS = "D"
        # NOT_INVOLVED = "X"
        # RESOLVED = "Z"
        # NO_CONFLICT = "N"

        active = qs.count()
        medical = qs.count(severity__gte=ConflictSeverity.MEDICALLY_SIGNIFICANT)
        waiting_for_triage = qs.filter(overall_status=DiscordanceReportNextStep.AWAITING_YOUR_TRIAGE).count()
        waiting_for_triage_medical = qs.filter(overall_status=DiscordanceReportNextStep.AWAITING_YOUR_TRIAGE, severity__gte=ConflictSeverity.MEDICALLY_SIGNIFICANT).count()
        waiting_for_amend = qs.filter(overall_status=DiscordanceReportNextStep.AWAITING_YOUR_AMEND).count()
        ready_to_discuss = qs.filter(overall_status=DiscordanceReportNextStep.TO_DISCUSS).count()

        return ConflictCategoriesUtils(
            active=active,
            medical=medical,
            waiting_for_triage=waiting_for_triage,
            waiting_for_triage_medical=waiting_for_triage_medical,
            waiting_for_amend=waiting_for_amend,
            ready_to_discuss=ready_to_discuss
        )
