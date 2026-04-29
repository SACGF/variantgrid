from django.core.management import BaseCommand

from annotation.models import ClinVarRecordCollection
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import TestingContextBucket
from classification.models import ClassificationGrouping, Overlap, OverlapContribution, ClassificationResultValue, \
    EvidenceKeyMap, EffectiveDateType, TriageStatus, EffectiveDate, TriageState
from classification.models.overlaps_enums import OverlapContributionStatus, OverlapEntrySourceTextChoices
from classification.services.overlaps_services import OverlapServices


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--migrate', required=False, action="store_true",
                            help="Migrates legacy data to new format - work in progress")
        parser.add_argument('--full_reset', required=False, action="store_true",
                            help="Deletes all Overlaps and OverlapContributions and creates them from scratch")


    def handle(self, *args, **options):
        if options["full_reset"]:
            self.full_reset(args, options)
        elif options["migrate"]:
            self.populate_status_change()

    def full_reset(self, *args, **options):
        Overlap.objects.all().delete()
        OverlapContribution.objects.all().delete()

        for cg in ClassificationGrouping.objects.iterator():
            OverlapServices.update_classification_grouping_overlap_contribution(cg)

        self.make_clinvar_expert_panel_contributions()

        print(f"Overlap Count = {Overlap.objects.count()}")
        print(f"Overlap Contribution Count = {OverlapContribution.objects.count()}")

    def populate_status_change(self):
        for overlap in Overlap.objects.iterator():
            latest_date = None
            for contribution in overlap.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING):
                if grouping := contribution.classification_grouping:
                    for mod in grouping.classification_modifications:
                        latest_mod_date = mod.created
                        if latest_date is None or latest_mod_date > latest_date:
                            latest_date = latest_mod_date
            if latest_date:
                overlap.overlap_status_change_timestamp = latest_date
                overlap.save(update_fields=["overlap_status_change_timestamp"])

    def make_clinvar_expert_panel_contributions(self):
        # only check already made ClinVarRecord collections in sync
        for clinvar_record_collection in ClinVarRecordCollection.objects.filter(
                max_stars__gte=CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE):
            if expert_panel := clinvar_record_collection.expert_panel:

                value = expert_panel.clinical_significance
                relevant_value = ClassificationResultValue.ONC_PATH and EvidenceKeyMap.clinical_significance_to_bucket().get(
                    value) is not None
                contribution_enum = OverlapContributionStatus.CONTRIBUTING if relevant_value else OverlapContributionStatus.NON_COMPARABLE_VALUE
                effective_date = EffectiveDate.from_datetime(expert_panel.date_last_evaluated or expert_panel.date_clinvar_updated, EffectiveDateType.CURATED)

                # FIXME this code is duplicated
                contribution, created = OverlapContribution.objects.update_or_create(
                    source=OverlapEntrySourceTextChoices.CLINVAR,
                    scv=expert_panel.record_id,
                    allele=clinvar_record_collection.allele,
                    classification_grouping=None,
                    value_type=ClassificationResultValue.ONC_PATH,
                    contribution_status=contribution_enum,
                    testing_context_bucket=TestingContextBucket.GERMLINE,
                    tumor_type_category=None,
                    defaults={
                        "value": value,
                        "effective_date": effective_date,
                    },
                    triage_state=TriageState(status=TriageStatus.NON_INTERACTIVE_THIRD_PARTY)
                )

                OverlapServices.link_overlap_contribution(contribution)
                for skew in contribution.overlapcontributionskew_set.select_related('overlap').all():
                    OverlapServices.recalc_overlap(skew.overlap)
