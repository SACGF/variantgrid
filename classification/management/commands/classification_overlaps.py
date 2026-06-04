from auditlog.context import disable_auditlog, set_extra_data
from django.core.management import BaseCommand
from django.db.models import F

from annotation.models import ClinVarRecordCollection, ClinVarRecord, EffectiveDate
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import TestingContextBucket, OverlapStatus, OverlapType
from classification.models import ClassificationGrouping, Overlap, OverlapContribution, ClassificationResultValue, \
    EvidenceKeyMap, EffectiveDateType, DiscordanceReport, DiscordanceReportTriageStatus, classification_flag_types, \
    ClassificationFlagTypes
from classification.enums.overlaps_enums import OverlapContributionStatus, OverlapEntrySourceTextChoices, TriageState, \
    TriageStatus
from classification.services.overlaps_services import OverlapServices
from snpdb.models import Allele


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--migrate', required=False, action="store_true",
                            help="Migrates legacy data to new format - work in progress")
        parser.add_argument('--full_reset', required=False, action="store_true",
                            help="Deletes all Overlaps and OverlapContributions and creates them from scratch")
        parser.add_argument("--recalc_skews", required=False, action="store_true", help="Updates what each lab's perspective of the Overlap is")
        parser.add_argument("--max_status", required=False, action="store_true",
                            help="Populates max status (not entirely accurate) but should distinguish between records that have had a discordance reports to the ones that haven't")


    def handle(self, *args, **options):
        if options["full_reset"]:
            self.full_reset(args, options)
        elif options["migrate"]:
            self.populate_overlap_change_date()
            self.migrate_discordance_reports(args, options)
        elif options["recalc_skews"]:
            self.recalc_skews()
        elif options["max_status"]:
            self.populate_max_status(args, options)
        else:
            print("Must choose full_reset or migrate")

    def recalc_skews(self):
        for overlap in Overlap.objects.all().iterator():
            OverlapServices.update_skews(overlap)

    def full_reset(self, *args, **options):
        Overlap.objects.all().delete()
        OverlapContribution.objects.all().delete()

        for cg in ClassificationGrouping.objects.iterator():
            OverlapServices.update_classification_grouping_overlap_contribution(cg, migration=True)

        self.make_clinvar_expert_panel_contributions()

        print(f"Overlap Count = {Overlap.objects.count()}")
        print(f"Overlap Contribution Count = {OverlapContribution.objects.count()}")

        self.populate_overlap_change_date()
        self.migrate_discordance_reports(args, options)
        self.populate_max_status()

    def populate_max_status(self, *args, **options):
        alleles_with_discordance_reports: dict[Allele, OverlapStatus] = {}
        for dr in DiscordanceReport.objects.filter(clinical_context__allele__isnull=False).select_related('clinical_context__allele').order_by('-created').iterator():
            allele = dr.clinical_context.allele
            if alleles_with_discordance_reports.get(allele) == OverlapStatus.MEDICALLY_SIGNIFICANT:
                continue

            if dr.is_medically_significant:
                alleles_with_discordance_reports[allele] = OverlapStatus.MEDICALLY_SIGNIFICANT
            elif dr.clinical_context.allele not in alleles_with_discordance_reports:
                alleles_with_discordance_reports[allele] = OverlapStatus.MAJOR_DIFFERENCES

        for allele, value in alleles_with_discordance_reports.items():
            Overlap.objects.filter(overlap_type=OverlapType.SINGLE_CONTEXT, value_type=ClassificationResultValue.ONC_PATH, allele=allele).update(overlap_max_ever_status=value)

        Overlap.objects.filter(overlap_max_ever_status__lt=F('overlap_status')).update(overlap_max_ever_status=F('overlap_status'))


    def migrate_discordance_reports(self, *args, **options):
        for discordance_report in DiscordanceReport.objects.all().order_by('created').iterator():
            allele = discordance_report.clinical_context.allele
            overlap = Overlap.objects.filter(
                allele=allele,
                testing_context_bucket=TestingContextBucket.GERMLINE,
            ).first()
            if not overlap:
                print(f"Could not find overlap for Discordance Report {discordance_report.pk}")
            else:
                # lab = models.ForeignKey(Lab, on_delete=CASCADE)
                #     triage_status = models.TextField(max_length=1, choices=DiscordanceReportTriageStatus.choices, default=DiscordanceReportTriageStatus.PENDING)
                #     note = models.TextField(null=True, blank=True)
                #     triage_date = models.DateField(null=True, blank=True)
                #     user = models.ForeignKey(User, null=True, blank=True, on_delete=PROTECT)

                for legacy_triage in discordance_report.discordancereporttriage_set.exclude(
                    triage_status=DiscordanceReportTriageStatus.PENDING,
                    note__isnull=True
                ):
                    overlap_contribution = OverlapContribution.objects.filter(
                        allele=allele,
                        testing_context_bucket=overlap.testing_context_bucket,
                        classification_grouping__lab=legacy_triage.lab).first()

                    triage_status: TriageStatus
                    try:
                        # there might be some legacy TriageStatus in the database
                        triage_status = TriageStatus(legacy_triage.triage_status)
                    except ValueError:
                        print(f"Found illegal triage status \"{legacy_triage.triage_status}\" in triage {legacy_triage.pk}")
                        continue

                    if triage_status == TriageStatus.PENDING and not legacy_triage.note:
                        # user did not update triage
                        continue

                    if not overlap_contribution:
                        print(f"Could not find OverlapContribution for DR_{discordance_report.pk} but there were modified triages against it")
                        continue

                    triage_state = TriageState(triage_status)
                    if triage_status == TriageStatus.REVIEWED_WILL_FIX:
                        found_pending_value = False
                        # see if we can get the flag from somewhere if it's outstanding
                        if target_grouping := overlap_contribution.classification_grouping:
                            if target_classification := target_grouping.latest_classification_modification:
                                if flag := target_classification.classification.flag_collection_safe.get_flag_of_type(classification_flag_types.classification_pending_changes):
                                    pending_value = (flag.data or {}).get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY)
                                    triage_state = TriageState(triage_status, pending_value)
                                    found_pending_value = True
                        if not found_pending_value:
                            print(f"Did not find pending change value for triage {legacy_triage.pk}, maybe it already changed")

                    has_change = False
                    if overlap_contribution.triage_state_obj != triage_state:
                        overlap_contribution.triage_state_obj = triage_state
                        has_change = True

                    if note := legacy_triage.note:
                        if note != overlap_contribution.comment_obj.text:
                            overlap_contribution.comment_obj = overlap_contribution.comment_obj.next_comment(note)
                            print(f"SETTING NOTE to {overlap_contribution.comment_obj}")
                            has_change = True

                    if has_change:
                        print(f"Updating triage cousin to DR_{discordance_report.pk} - on allele {overlap.allele:CA}")
                        with set_extra_data({
                            "actor": legacy_triage.user,
                            "remote_addr": None,
                            "remote_port": None,
                            "timestamp": legacy_triage.modified,
                            "migration": True
                        }):
                            overlap_contribution.save()
                    else:
                        print(f"Triage already up to date - possibly multiple discordances for same record")


    def populate_overlap_change_date(self):
        # timestamp on overlaps
        with disable_auditlog():
            for overlap in Overlap.objects.all().iterator():
                # note we're looking for the latest published date of a classification here
                # as the upload date is when a discordance would occur
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

            # dates on overlap contributions
            for overlap_contribution in OverlapContribution.objects.filter(effective_date__date=None).iterator():
                if grouping := overlap_contribution.classification_grouping:
                    if latest_modification := grouping.latest_classification_modification:
                        date_check = latest_modification.curated_date_check
                        overlap_contribution.effective_date_obj = date_check.to_effective_date
                        print(f"Setting effective date to {date_check.to_effective_date}")
                        print(overlap_contribution.effective_date)
                        overlap_contribution.save(update_fields=["effective_date"])
                elif scv := overlap_contribution.scv:
                    if clinvar_record := ClinVarRecord.objects.filter(record_id=scv).first():
                        overlap_contribution.effective_date_obj = clinvar_record.effective_date
                        print(f"Setting effective date to {date_check.to_effective_date}")
                        print(overlap_contribution.effective_date)
                        overlap_contribution.save(update_fields=["effective_date"])

    def make_clinvar_expert_panel_contributions(self):
        # only check already made ClinVarRecord collections in sync
        for clinvar_record_collection in ClinVarRecordCollection.objects.filter(
                max_stars__gte=CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE, allele__isnull=False):
            expert_panel: ClinVarRecord
            if expert_panel := clinvar_record_collection.expert_panel:

                value = expert_panel.clinical_significance
                relevant_value = ClassificationResultValue.ONC_PATH and EvidenceKeyMap.clinical_significance_to_bucket().get(
                    value) is not None
                contribution_enum = OverlapContributionStatus.CONTRIBUTING if relevant_value else OverlapContributionStatus.NON_COMPARABLE_VALUE
                effective_date = EffectiveDate.from_datetime(expert_panel.date_last_evaluated or expert_panel.date_clinvar_updated, EffectiveDateType.CURATED)

                # WARNING: Do ClinVarRecords get wiped and replaced? but even if so, we will keep this created (in the migration) as the first
                # date
                with set_extra_data({"timestamp": expert_panel.created, "migration": True}):
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
                            "effective_date": effective_date.to_dict(),
                        },
                        triage_state=TriageState(status=TriageStatus.NON_INTERACTIVE_THIRD_PARTY).to_dict()
                    )

                    OverlapServices.link_overlap_contribution(contribution)
                    for skew in contribution.overlapcontributionskew_set.select_related('overlap').all():
                        OverlapServices.recalc_overlap(skew.overlap)
