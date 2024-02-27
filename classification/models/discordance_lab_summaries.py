from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from classification.enums import SpecialEKeys
from classification.models import ClassificationLabSummary, DiscordanceReport, ClassificationLabSummaryEntry, \
    classification_flag_types, ClassificationFlagTypes, DiscordanceReportClassification, DiscordanceReportTriage
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab


@dataclass(frozen=True)
class DiscordanceLabSummary(ClassificationLabSummary):
    drcs: list[DiscordanceReportClassification]
    triage: Optional[DiscordanceReportTriage] = None

    @cached_property
    def embedded(self):
        if triage := self.triage:
            from classification.views.discordance_report_triage_view import DiscordanceReportTriageView
            return DiscordanceReportTriageView.lazy_render(triage)

    def _with_triage(self, triage: Optional[DiscordanceReportTriage]) -> 'DiscordanceLabSummary':
        return DiscordanceLabSummary(
            group=self.group,
            is_internal=self.is_internal,
            count=self.count,
            drcs=self.drcs,
            triage=triage
        )

    @staticmethod
    def for_discordance_report(discordance_report: DiscordanceReport, perspective: LabPickerData) -> list['DiscordanceLabSummary']:
        group_counts: dict[ClassificationLabSummaryEntry, list[DiscordanceReportClassification]] = defaultdict(list)
        for drc in DiscordanceReportClassification.objects.filter(report=discordance_report).select_related(
                'classification_original',
                'classification_original__classification',
                'classification_original__classification__lab',
                'classification_original__classification__lab__organization',
                'classification_final'
        ):
            clinical_significance_from = drc.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            clinical_significance_to: Optional[str] = None
            pending: bool = False
            if drc.withdrawn_effective:
                clinical_significance_to = 'withdrawn'
            else:
                if discordance_report.is_latest:
                    # classification is still outstanding, check to see if there pending changes
                    if flag := drc.classification_original.classification.flag_collection.get_open_flag_of_type(
                            classification_flag_types.classification_pending_changes):
                        clinical_significance_to = flag.data.get(
                            ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY) or "unknown"
                        pending = True

                if not clinical_significance_to:
                    clinical_significance_to = drc.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

            group_counts[ClassificationLabSummaryEntry(
                lab=drc.classification_original.classification.lab,
                clinical_significance_from=clinical_significance_from,
                clinical_significance_to=clinical_significance_to,
                somatic_clinical_significance=drc.classification_original.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE),
                amp_level=drc.classification_original.amp_level,
                pending=pending
            )].append(drc)

        dlses: list[DiscordanceLabSummary] = list(sorted([DiscordanceLabSummary(
            group=group,
            is_internal=group.lab in perspective.labs_if_not_admin,
            count=len(drcs),
            drcs=drcs
        ) for group, drcs in group_counts.items()]))

        triage_by_lab: dict[Lab, DiscordanceReportTriage] = {}
        for triage in discordance_report.discordancereporttriage_set.select_related('lab').all():
            triage_by_lab[triage.lab] = triage

        # apply triages to the first (non withdrawn row) for each lab
        # in case a lab has both withdrawn and non-withdrawn records
        with_triages = []
        for dl in dlses:
            if not dl.clinical_significance_to == 'withdrawn':
                dl = dl._with_triage(triage_by_lab.pop(dl.lab, None))
            with_triages.append(dl)

        if triage_by_lab:
            # there are some triages for purely withdrawn labs
            for dl in with_triages:
                if not dl.triage:
                    dl = dl._with_triage(triage_by_lab.pop(dl.lab, None))
            with_triages.append(dl)

        return with_triages
