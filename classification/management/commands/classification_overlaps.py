from datetime import datetime
from typing import Optional

from auditlog.context import disable_auditlog, set_extra_data
from auditlog.models import LogEntry
from django.contrib.contenttypes.models import ContentType
from django.core.management import BaseCommand
from django.db import transaction
from django.db.models import F, OuterRef, Q, Subquery
from django.db.models.functions import Coalesce

from annotation.models import ClinVarRecordCollection
from annotation.models.data_enums import EffectiveDate
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import TestingContextBucket, OverlapStatus, OverlapType, SpecialEKeys, \
    OverlapEntrySourceTextChoices, ShareLevel, AlleleOriginBucket
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import ClassificationGrouping, Overlap, OverlapContribution, ClassificationResultValue, \
    DiscordanceReport, DiscordanceReportTriageStatus, classification_flag_types, \
    ClassificationFlagTypes, ClassificationModification, Classification, ClassificationGroupingEntry, \
    OverlapContributionNextStep
from classification.enums.overlaps_enums import OverlapContributionStatus, TriageState, \
    TriageStatus
from classification.services.overlap_calculator import overlap_calculator_for_value_type
from classification.services.overlaps_services import OverlapServices, OverlapPageDetails
from flags.models import Flag, FlagStatus
from library.django_utils import chunked_queryset
from snpdb.models import Allele


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--flags', required=False, action='store_true'),
        parser.add_argument('--cache', required=False, action='store_true'),
        parser.add_argument('--dates', required=False, action='store_true'),
        parser.add_argument('--recalc', required=False, action="store_true"),
        parser.add_argument('--full_reset', required=False, action="store_true",
                            help="Deletes all Overlaps and OverlapContributions and creates them from scratch")
        parser.add_argument("--recalc_skews", required=False, action="store_true", help="Updates what each lab's perspective of the Overlap is")
        parser.add_argument("--max_status", required=False, action="store_true",
                            help="Populates max status (not entirely accurate) but should distinguish between records that have had a discordance reports to the ones that haven't")

    def handle(self, *args, **options):
        if options["full_reset"]:
            # rebuilding from scratch makes every existing discordance look new
            with OverlapServices.discordance_notifications_suppressed():
                self.full_reset(args, options)
            return

        if options['flags']:
            self.check_pending_change_flags()
            self.recalc_overlaps()

        if options['recalc']:
            self.recalc_overlaps()

        if options['cache']:
            bulk_update = []
            for overlap in Overlap.objects.all().iterator():
                overlap.cached_overlap_state_obj = overlap.derived_overlap_state
                bulk_update.append(overlap)
            print(f"Bulk updating {len(bulk_update)}")
            update_count = Overlap.objects.bulk_update(bulk_update, fields=["cached_overlap_state"])
            print(f"Updated {update_count}")

        if options["recalc_skews"]:
            self.recalc_skews()

        if options["max_status"]:
            self.populate_max_status(args, options)

        if options["dates"]:
            self.populate_overlap_change_date()

    def recalc_overlaps(self):
        print("Recalcing Overlaps")
        OverlapServices.recalc_overlaps_batch(Overlap.objects.order_by("pk").values_list("pk", flat=True))

    def recalc_skews(self):
        print("Recalcing just skews")
        for overlap in Overlap.objects.all().iterator():
            OverlapServices.update_next_steps(overlap)

    def full_reset(self, *args, **options):
        print("Full Reset")
        with disable_auditlog():
            Overlap.objects.all().delete()
            OverlapContribution.objects.all().delete()

        self.populate_overlap_history()
        self.migrate_discordance_reports(args, options)
        self.populate_max_status()
        self.check_pending_change_flags()
        self.make_clinvar_expert_panel_contributions()
        self.recalc_overlaps()
        self.populate_overlap_change_date()

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

    def check_pending_change_flags(self):
        for flag in Flag.objects.filter(flag_type=classification_flag_types.classification_pending_changes):
            if flag.resolution.status == FlagStatus.OPEN:
                if pending_value := (flag.data or {}).get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY):
                    if classification := Classification.objects.filter(flag_collection=flag.collection).first():
                        if grouping := ClassificationGroupingEntry.grouping_for(classification):
                            # TODO see if we can use audit log to time this as part of the original
                            triage: OverlapContribution = grouping.contribution_for(ClassificationResultValue.ONC_PATH)
                            new_triage_state = TriageState(TriageStatus.REVIEWED_WILL_FIX, pending_value)
                            if triage.triage_state_obj != new_triage_state:
                                print("Updating triage with pending value")
                                with disable_auditlog():
                                    triage.triage_state_obj = new_triage_state
                                    triage.save()


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
        print("Recalcing Dates")
        # an overlap last changed with the latest audit log entry of any of its contributions
        contribution_ids = OverlapContributionNextStep.objects.filter(overlap=OuterRef(OuterRef("pk"))).values("contribution_id")
        latest_log_entry = LogEntry.objects.filter(
            content_type=ContentType.objects.get_for_model(OverlapContribution),
            object_id__in=contribution_ids,
        ).order_by("-timestamp").values("timestamp")[:1]
        Overlap.objects.update(overlap_status_change_timestamp=Coalesce(Subquery(latest_log_entry), F("overlap_status_change_timestamp")))

        # or, for germline, with a discordance report on its allele (or a review of one)
        allele_discordance_dates: dict[int, datetime] = {}
        for dr in DiscordanceReport.objects.filter(clinical_context__allele__isnull=False).select_related("clinical_context"):
            allele_id = dr.clinical_context.allele_id
            dates = [dr.modified] + [review.modified for review in dr.reviews_all()]
            if previous_date := allele_discordance_dates.get(allele_id):
                dates.append(previous_date)
            allele_discordance_dates[allele_id] = max(dates)

        for allele_id, discordance_date in allele_discordance_dates.items():
            Overlap.objects.filter(
                Q(overlap_status_change_timestamp__isnull=True) | Q(overlap_status_change_timestamp__lt=discordance_date),
                allele_id=allele_id,
                testing_context_bucket=TestingContextBucket.GERMLINE,
            ).update(overlap_status_change_timestamp=discordance_date)

    def populate_overlap_history(self):
        grouping_qs = ClassificationGrouping.objects.select_related(
            "allele_origin_grouping__allele", "lab__organization", "latest_classification_modification")
        for grouping_chunk_qs in chunked_queryset(grouping_qs, 500):
            with transaction.atomic():
                for grouping in grouping_chunk_qs:
                    self._populate_grouping_overlap_history(grouping)

    @staticmethod
    def _populate_grouping_overlap_history(grouping: ClassificationGrouping):
        value_types = [ClassificationResultValue.ONC_PATH]
        if grouping.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
            value_types = [ClassificationResultValue.ONC_PATH, ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE]
        for value_type in value_types:
            calc = overlap_calculator_for_value_type(value_type)
            all_classifications = grouping.classificationgroupingentry_set.values_list("classification", flat=True)
            all_modifications = ClassificationModification.objects.filter(classification_id__in=all_classifications, published=True, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).order_by('created')

            overlap_contribution: Optional[OverlapContribution] = None
            last_value: Optional[str] = None
            last_curation_date: Optional[EffectiveDate] = None
            for modification in all_modifications.iterator():
                this_value = None
                if value_type == ClassificationResultValue.ONC_PATH:
                    this_value = modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                elif value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
                    this_value = modification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)

                if this_value is None:
                    continue # don't log uncurated variants

                this_curation_date = modification.curated_date_check.to_effective_date
                if last_curation_date and last_curation_date > this_curation_date:
                    # don't go backwards with curation dates
                    continue

                this_modification_date = modification.created

                if this_value != last_value or this_curation_date != last_curation_date:
                    contribution_status = OverlapContributionStatus.CONTRIBUTING if calc.is_comparable_value(this_value) else OverlapContributionStatus.NON_COMPARABLE_VALUE
                    audit_context = {"migration": True, "timestamp": this_modification_date}
                    with set_extra_data(audit_context):
                        if not overlap_contribution:
                            overlap_contribution = OverlapContribution.objects.create(
                                source=OverlapEntrySourceTextChoices.CLASSIFICATION,
                                allele=grouping.allele,
                                classification_grouping=grouping,
                                testing_context_bucket=grouping.testing_context,
                                tumor_type_category=grouping.tumor_type_category,
                                value_type=value_type,
                                value=this_value,
                                contribution_status=contribution_status,
                                effective_date=this_curation_date.to_dict()
                            )
                            OverlapServices._link_overlap_contribution(overlap_contribution)
                        else:
                            overlap_contribution.value = this_value
                            overlap_contribution.contribution_status = contribution_status
                            overlap_contribution.effective_date_obj = this_curation_date
                            overlap_contribution.save()

                    last_value = this_value
                    last_curation_date = this_curation_date

        # This will make it so we still create OverlapContributions for unshared records
        # as well as making sure the data we end up with matches the migration script
        OverlapServices.update_classification_grouping_overlap_contribution(grouping, migration=True, recalc_overlaps=False)

    def make_clinvar_expert_panel_contributions(self):
        # only check already made ClinVarRecord collections in sync
        for clinvar_record_collection in ClinVarRecordCollection.objects.filter(
                max_stars__gte=CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE, allele__isnull=False):
            if clinvar_record_collection.expert_panel is not None:
                OverlapServices.update_clinvar_overlap_contribution(clinvar_record_collection, migration=True, recalc_overlaps=False)
