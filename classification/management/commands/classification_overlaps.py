from auditlog.context import disable_auditlog, set_extra_data
from django.core.management import BaseCommand
from django.db.models import F

from annotation.models import ClinVarRecordCollection
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import TestingContextBucket, OverlapStatus, OverlapType
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import ClassificationGrouping, Overlap, OverlapContribution, ClassificationResultValue, \
    DiscordanceReport, DiscordanceReportTriageStatus, classification_flag_types, \
    ClassificationFlagTypes
from classification.enums.overlaps_enums import OverlapContributionStatus, TriageState, \
    TriageStatus
from classification.services.overlaps_services import OverlapServices
from snpdb.models import Allele


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--recalc', required=False, action="store_true"),
        parser.add_argument('--migrate', required=False, action="store_true",
                            help="Migrates legacy data to new format - work in progress")
        parser.add_argument('--full_reset', required=False, action="store_true",
                            help="Deletes all Overlaps and OverlapContributions and creates them from scratch")
        parser.add_argument("--recalc_skews", required=False, action="store_true", help="Updates what each lab's perspective of the Overlap is")
        parser.add_argument("--max_status", required=False, action="store_true",
                            help="Populates max status (not entirely accurate) but should distinguish between records that have had a discordance reports to the ones that haven't")

    def handle(self, *args, **options):
        if options['recalc']:
            self.recalc_overlaps()
        elif options["full_reset"]:
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

    def recalc_overlaps(self):
        for overlap in Overlap.objects.all().iterator():
            OverlapServices.recalc_overlap(overlap)
            OverlapServices.update_skews(overlap)

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

                is_continued_discordance = discordance_report.is_latest and discordance_report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE
                if is_continued_discordance:
                    print(f"*** CONTINUED DISCORDANCE for OV_{overlap.pk} ***")

                legacy_triage_qs = discordance_report.discordancereporttriage_set.all()
                if not is_continued_discordance:
                    legacy_triage_qs = legacy_triage_qs.exclude(triage_status=DiscordanceReportTriageStatus.PENDING)

                for legacy_triage in legacy_triage_qs:
                    has_stuff = legacy_triage.triage_status != DiscordanceReportTriageStatus.PENDING

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

                    if has_stuff and not overlap_contribution:
                        print(f"Could not find OverlapContribution for DR_{discordance_report.pk} but there were modified triages against it")
                        continue

                    if triage_status == TriageStatus.PENDING and not legacy_triage.note and not is_continued_discordance:
                        # user did not update triage, and this isn't a continued discordance so nothing to report
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

                    if is_continued_discordance:
                        print(f"Migrating continued discordance")
                        overlap_contribution.review_agreed_value = overlap_contribution.effective_value
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
                if is_continued_discordance:
                    OverlapServices.recalc_overlap(overlap)

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
                        print(overlap_contribution.effective_date)
                        overlap_contribution.save(update_fields=["effective_date"])

    def make_clinvar_expert_panel_contributions(self):
        # only check already made ClinVarRecord collections in sync
        for clinvar_record_collection in ClinVarRecordCollection.objects.filter(
                max_stars__gte=CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE, allele__isnull=False):
            if clinvar_record_collection.expert_panel is not None:
                OverlapServices.update_clinvar_overlap_contribution(clinvar_record_collection, migrate=True)
