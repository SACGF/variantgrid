from django.core.management import BaseCommand

from annotation.models import ClinVarRecordCollection
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import TestingContextBucket
from classification.models import ClassificationGrouping, Overlap, OverlapContribution, ClassificationResultValue, \
    EvidenceKeyMap
from classification.models.overlaps_enums import OverlapContributionStatus, OverlapEntrySourceTextChoices
from classification.services.overlaps_services import OverlapServices


class Command(BaseCommand):

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        Overlap.objects.all().delete()
        OverlapContribution.objects.all().delete()

        for cg in ClassificationGrouping.objects.iterator():
            OverlapServices.update_classification_grouping_overlap_contribution(cg)

        self.make_clinvar_expert_panel_contributions()

        print(f"Overlap Count = {Overlap.objects.count()}")
        print(f"Overlap Contribution Count = {OverlapContribution.objects.count()}")

    def make_clinvar_expert_panel_contributions(self):
        # only check already made ClinVarRecord collections in sync
        for clinvar_record_collection in ClinVarRecordCollection.objects.filter(
                max_stars__gte=CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE):
            if expert_panel := clinvar_record_collection.expert_panel:
                value = expert_panel.clinical_significance
                relevant_value = ClassificationResultValue.ONC_PATH and EvidenceKeyMap.clinical_significance_to_bucket().get(
                    value) is not None
                contribution_enum = OverlapContributionStatus.CONTRIBUTING if relevant_value else OverlapContributionStatus.NON_COMPARABLE_VALUE
                effective_date = expert_panel.date_last_evaluated or expert_panel.date_clinvar_updated

                contribution, _ = OverlapContribution.objects.get_or_create(
                    source=OverlapEntrySourceTextChoices.CLINVAR,
                    scv=expert_panel.record_id,
                    allele=clinvar_record_collection.allele,
                    classification_grouping=None,
                    value_type=ClassificationResultValue.ONC_PATH,
                    value=value,
                    contribution=contribution_enum,
                    testing_context_bucket=TestingContextBucket.GERMLINE,
                    tumor_type_category=None,
                    effective_date=effective_date
                )
                OverlapServices.link_overlap_contribution(contribution)
                for skew in contribution.overlapcontributionskew_set.select_related('overlap').all():
                    OverlapServices.recalc_overlap(skew.overlap)
