from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, Dict, List

from classification.enums import SpecialEKeys
from classification.models import ClassificationLabSummary, DiscordanceReport, ClassificationLabSummaryEntry, \
    classification_flag_types, ClassificationFlagTypes, DiscordanceReportClassification
from snpdb.lab_picker import LabPickerData


@dataclass(frozen=True)
class DiscordanceLabSummary(ClassificationLabSummary):
    drcs: List[DiscordanceReportClassification]

    @staticmethod
    def for_discordance_report(discordance_report: DiscordanceReport, perspective: LabPickerData) -> List['DiscordanceLabSummary']:
        group_counts: Dict[ClassificationLabSummaryEntry, List[DiscordanceReportClassification]] = defaultdict(list)
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
                pending=pending
            )].append(drc)

        return sorted([DiscordanceLabSummary(
            group=group,
            is_internal=group.lab in perspective.labs_if_not_admin,
            count=len(drcs),
            drcs=drcs
        ) for group, drcs in group_counts.items()])