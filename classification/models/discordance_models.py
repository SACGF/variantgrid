import typing
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from enum import Enum
from typing import Set, Optional, List, Dict, Tuple, Any, Iterable

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import PROTECT, CASCADE
from django.urls.base import reverse
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from lazy import lazy
from more_itertools import first

from classification.enums.classification_enums import SpecialEKeys, ClinicalSignificance
from classification.enums.discordance_enums import DiscordanceReportResolution, ContinuedDiscordanceReason
from classification.models.classification import ClassificationModification, Classification
from classification.models.classification_lab_summaries import ClassificationLabSummaryEntry, ClassificationLabSummary
from classification.models.clinical_context_models import ClinicalContext
from classification.models.flag_types import classification_flag_types, ClassificationFlagTypes
from flags.models.enums import FlagStatus
from flags.models.models import FlagComment
from genes.hgvs import CHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab, GenomeBuild

discordance_change_signal = django.dispatch.Signal()  # args: "discordance_report", "cause"


class NotifyLevel(str, Enum):
    NEVER_NOTIFY = "never-notify"
    NOTIFY_IF_CHANGE = "notify-if-change"
    ALWAYS_NOTIFY = "always-notify"


class DiscordanceReport(TimeStampedModel):

    resolution = models.TextField(default=DiscordanceReportResolution.ONGOING, choices=DiscordanceReportResolution.CHOICES, max_length=1, null=True, blank=True)
    # TODO remove continued discordane reason, it should be redundant to notes
    continued_discordance_reason = models.TextField(choices=ContinuedDiscordanceReason.CHOICES, max_length=1, null=True, blank=True)
    notes = models.TextField(null=False, blank=True, default='')

    clinical_context = models.ForeignKey(ClinicalContext, on_delete=PROTECT)

    # report started date is a bit redundant compared to TimeStampModel's created
    # but as it has more of a business implication than a technical one I felt it
    # should be its own field
    report_closed_by = models.ForeignKey(User, on_delete=PROTECT, null=True)
    report_started_date = models.DateTimeField(auto_now_add=True)
    report_completed_date = models.DateTimeField(null=True, blank=True)

    cause_text = models.TextField(null=False, blank=True, default='')
    resolved_text = models.TextField(null=False, blank=True, default='')

    class LabInvolvement(int, Enum):
        WITHDRAWN = 1
        ACTIVE = 2

    def get_absolute_url(self):
        return reverse('discordance_report', kwargs={"discordance_report_id": self.pk})

    def __str__(self):
        return f"({self.pk}) Discordance Report ({self.clinical_context.allele:CA}) {self.get_resolution_display() or 'Discordant'}"

    @property
    def metrics_logging_key(self) -> Tuple[str, Any]:
        return "discordance_report_id", self.pk

    @property
    def status(self) -> str:
        if self.resolution:
            return self.get_resolution_display()
        return "Active Discordance"

    @property
    def is_important(self):
        return self.resolution is None or self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE

    @property
    def is_active(self) -> bool:
        return self.resolution is None

    @property
    def resolution_text(self) -> str:
        if self.resolution == DiscordanceReportResolution.CONCORDANT:
            return 'Concordant'
        if self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
            return 'Continued Discordance'
        if self.is_pending_concordance:
            return 'Pending Concordance'
        if self.resolution == DiscordanceReportResolution.ONGOING:
            return 'In Discordance'
        return ''

    @transaction.atomic
    def close(self, expected_resolution: Optional[str] = None, cause_text: str = ''):
        """
        @param expected_resolution will be calculated, but if provided a ValueError will be raised if it doesn't match
        @param cause_text reason this discordance is being closed
        """
        if self.clinical_context.discordance_status.is_discordant:
            self.resolution = DiscordanceReportResolution.CONTINUED_DISCORDANCE
        else:
            self.resolution = DiscordanceReportResolution.CONCORDANT

        if expected_resolution and expected_resolution != self.resolution:
            raise ValueError(f'Expected to close Discordance Report as {expected_resolution} but was {self.resolution}')

        self.resolved_text = cause_text
        self.report_completed_date = datetime.now()
        self.save()

        for drc in DiscordanceReportClassification.objects.filter(report=self):  # type: DiscordanceReportClassification
            drc.close()

        discordance_change_signal.send(DiscordanceReport, discordance_report=self, cause=cause_text)

    @transaction.atomic
    def update(self, cause_text: str = '', notify_level: NotifyLevel = NotifyLevel.NOTIFY_IF_CHANGE):
        if not self.is_active:
            raise ValueError('Cannot update non active Discordance Report')

        if not self.id:
            self.save()

        existing_vms = set()
        existing_labs = set()
        for drc in DiscordanceReportClassification.objects.filter(report=self):
            existing_vms.add(drc.classification_original.classification_id)
            existing_labs.add(drc.classification_original.classification.lab)

        newly_added_labs: Set[Lab] = set()
        for vcm_id in self.clinical_context.classifications_qs.values_list('id', flat=True):
            if vcm_id in existing_vms:
                existing_vms.remove(vcm_id)
            else:
                vcm = ClassificationModification.objects.get(is_last_published=True, classification=vcm_id)
                if vcm.classification.lab not in existing_labs:
                    newly_added_labs.add(vcm.classification.lab)

                DiscordanceReportClassification(
                    report=self,
                    classification_original=vcm
                ).save()

        # contents of existing_vms are now presumably records that changed clinical groupings

        if not self.clinical_context.discordance_status.is_discordant:
            # hey we've reached concordance, let's close this
            # close will fire off a notification event
            self.close(expected_resolution=DiscordanceReportResolution.CONCORDANT, cause_text=cause_text)
        else:
            if notify_level == NotifyLevel.NEVER_NOTIFY:
                pass
            elif notify_level == NotifyLevel.ALWAYS_NOTIFY:
                discordance_change_signal.send(DiscordanceReport, discordance_report=self, cause=cause_text)
            elif newly_added_labs:  # change is significant
                newly_added_labs_str = ", ".join(str(lab) for lab in newly_added_labs)
                discordance_change_signal.send(DiscordanceReport, discordance_report=self, cause=f"{cause_text} and newly added labs {newly_added_labs_str}")

    @property
    def all_actively_involved_labs(self) -> Set[Lab]:
        return {lab for lab, status in self.involved_labs.items() if status == DiscordanceReport.LabInvolvement.ACTIVE}

    @lazy
    def involved_labs(self) -> Dict[Lab, LabInvolvement]:
        lab_status: Dict[Lab, DiscordanceReport.LabInvolvement] = dict()
        for drc in DiscordanceReportClassification.objects.filter(report=self):
            effective_c: Classification = drc.classification_effective.classification
            lab = effective_c.lab
            status = DiscordanceReport.LabInvolvement.WITHDRAWN if effective_c.withdrawn else DiscordanceReport.LabInvolvement.ACTIVE
            lab_status[lab] = max(lab_status.get(lab, 0), status)
        return lab_status

    def user_is_involved(self, user: User) -> bool:
        if user.is_superuser:
            return True
        user_labs = set(Lab.valid_labs_qs(user, admin_check=True))
        return bool(user_labs.intersection(self.involved_labs.keys()))

    @transaction.atomic
    def reopen_continued_discordance(self, cause: str = '') -> 'DiscordanceReport':
        if not self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
            raise ValueError(f'Should only call create_new_report from the latest report, only if resolution = {DiscordanceReportResolution.CONTINUED_DISCORDANCE}')

        report = DiscordanceReport(clinical_context=self.clinical_context, cause_text=cause)
        report.save()
        # there used to be code that would create a new discordance report at the state the old "continued discordance" report was in
        # but that's very confusing as the data won't match the created date aka "discordance detected date"
        report.update(cause_text=cause, notify_level=NotifyLevel.ALWAYS_NOTIFY)
        return report

    @property
    def has_significance_changed(self):
        existing_vms = {}
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            if dr.clinical_context_final == self.clinical_context:
                existing_vms[dr.classification_final.classification.id] = dr.classification_final.get(SpecialEKeys.CLINICAL_SIGNIFICANCE, None)

        for vcm in self.clinical_context.classification_modifications:
            vc_id = vcm.classification.id
            if vc_id not in existing_vms:
                return True  # new classification, previous state is not relevant
            if existing_vms[vc_id] != vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                return True  # clinical significance has changed
            existing_vms.pop(vc_id, None)

        if existing_vms:
            return True
        return False

    @property
    def should_reopen_continued_discordance(self):
        if not self.clinical_context.discordance_status.is_discordant:
            # what was once "continued discordance" is concordant, re-open so we can instantly close it
            return True

        existing_labs: Set[Lab] = set()
        for drc in self.discordancereportclassification_set.all():
            existing_labs.add(drc.classification_original.classification.lab)

        # a new lab has submitted since continued discordance
        for cc in self.clinical_context.classifications_qs:
            if cc.lab not in existing_labs:
                return True
        return False

    @staticmethod
    def latest_report(clinical_context: ClinicalContext) -> 'DiscordanceReport':
        return DiscordanceReport.objects.filter(clinical_context=clinical_context).order_by('-created').first()

    @property
    def is_latest(self):
        return DiscordanceReport.latest_report(self.clinical_context) == self

    @staticmethod
    def update_latest(clinical_context: ClinicalContext, cause: str, update_flags: bool):
        try:
            latest_report = DiscordanceReport.latest_report(clinical_context=clinical_context)
            if latest_report:
                if latest_report.is_active:
                    latest_report.update(cause_text=cause)
                    return latest_report
                if latest_report.resolution == DiscordanceReportResolution.CONCORDANT:
                    # most recent report ended with discordance, can start fresh
                    pass
                elif latest_report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
                    # create a new report if data has changed significantly
                    if latest_report.should_reopen_continued_discordance:
                        return latest_report.reopen_continued_discordance(cause=f'Change after previous report marked as continued discordance - {cause}')
                    return None

            if clinical_context.is_discordant():
                report = DiscordanceReport(clinical_context=clinical_context, cause_text=cause)
                report.update(cause_text=cause)
                return report

            # no report has been made
            return None
        finally:
            if update_flags:
                DiscordanceReport.apply_flags_to_context(clinical_context)

    @staticmethod
    def apply_flags_to_context(clinical_context: ClinicalContext):
        # this is now the one function that applies and closes discordance flags
        # to classifications and clinical contexts
        discordant_classifications: Set[int] = set()
        is_in_discordance = False
        if latest_report := DiscordanceReport.latest_report(clinical_context):
            if latest_report.is_important:
                discordant_classifications = latest_report.actively_discordant_classification_ids()
                is_in_discordance = True

        # apply flags to clinical context
        if is_in_discordance:
            clinical_context.flag_collection_safe.get_or_create_open_flag_of_type(
                flag_type=classification_flag_types.clinical_context_discordance
            )
        else:
            clinical_context.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.clinical_context_discordance
            )

        # apply flags to classifications
        for vc in Classification.objects.filter(clinical_context=clinical_context):
            discordant = vc.id in discordant_classifications
            if discordant:
                vc.flag_collection_safe.get_or_create_open_flag_of_type(
                    flag_type=classification_flag_types.discordant
                )
            else:
                vc.flag_collection_safe.close_open_flags_of_type(
                    flag_type=classification_flag_types.discordant
                )

    def apply_flags(self):
        if self != DiscordanceReport.latest_report(self.clinical_context):
            raise ValueError("Can only apply_flags on the latest discordance report for any clinical context")
        DiscordanceReport.apply_flags_to_context(self.clinical_context)

    def actively_discordant_classification_ids(self) -> Set[int]:
        if not self.is_important:
            return set()
        classifications = set()
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            if dr.clinical_context_effective == self.clinical_context and not dr.withdrawn_effective:
                classifications.add(dr.classification_original.classification.id)
        return classifications

    @lazy
    def all_classification_modifications(self) -> List[ClassificationModification]:
        vcms: List[ClassificationModification] = list()
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            vcms.append(dr.classification_effective)
        return vcms

    def all_c_hgvs(self, genome_build: Optional[GenomeBuild] = None) -> List[CHGVS]:
        if not genome_build:
            genome_build = GenomeBuildManager.get_current_genome_build()
        c_hgvs = set()
        for cm in self.all_classification_modifications:
            c_hgvs.add(cm.c_hgvs_best(genome_build))
        return sorted(c_hgvs)

    @property
    def is_pending_concordance(self):
        return self.clinical_context.discordance_status.pending_concordance


# TODO all the below classes are utilites, consider moving them out


class DiscordanceReportRowData:

    # TODO merge this with DiscordanceReportViewTemplate?

    def __init__(self, discordance_report: DiscordanceReport, perspective: LabPickerData):
        self.discordance_report = discordance_report
        self.perspective = perspective

    @lazy
    def all_actively_involved_labs(self):
        return self.discordance_report.all_actively_involved_labs

    def all_actively_involved_labs_ids(self) -> str:
        return ";".join(str(lab.pk) for lab in self.all_actively_involved_labs)

    @property
    def is_requiring_attention(self) -> bool:
        return (self.perspective.is_admin_mode or bool(self.all_actively_involved_labs.intersection(self.perspective.your_labs))) \
               and self.discordance_report.resolution is None and not self.discordance_report.is_pending_concordance

    @property
    def is_valid_including_withdraws(self) -> bool:
        return self.perspective.is_admin_mode or bool(set(self.discordance_report.involved_labs.keys()).intersection(self.perspective.your_labs))

    @property
    def other_labs(self) -> Set[Lab]:
        return self.all_actively_involved_labs - self.perspective.labs_if_not_admin

    @property
    def is_internal(self):
        return not self.perspective.is_admin_mode and not self.other_labs

    @property
    def _cm_candidate(self) -> ClassificationModification:
        return first(self._cm_candidates)

    @property
    def _cm_candidates(self) -> Iterable[ClassificationModification]:
        yielded_any = False
        for cm_candidate in self.discordance_report.all_classification_modifications:
            if cm_candidate.classification.lab in self.perspective.your_labs:
                yielded_any = True
                yield cm_candidate
        if not yielded_any:
            yield self.discordance_report.all_classification_modifications[0]

    @property
    def id(self):
        return self.discordance_report.id

    @property
    def date_detected(self) -> datetime:
        return self.discordance_report.report_started_date

    @property
    def date_detected_str(self) -> str:
        date_str = f"{self.date_detected:%Y-%m-%d}"
        if (now() - self.date_detected) <= timedelta(days=1):
            date_str = f"{date_str} (NEW)"
        return date_str

    @property
    def c_hgvs(self) -> CHGVS:
        return self._cm_candidate.c_hgvs_best(genome_build=self.perspective.genome_build)

    @property
    def c_hgvses(self) -> List[CHGVS]:
        return sorted({candidate.c_hgvs_best(genome_build=self.perspective.genome_build) for candidate in self._cm_candidates})

    @property
    def is_pending_concordance(self):
        return self.discordance_report.is_pending_concordance

    @lazy
    def lab_significances(self) -> List[ClassificationLabSummary]:
        group_counts: Dict[ClassificationLabSummaryEntry, int] = defaultdict(int)
        for drc in DiscordanceReportClassification.objects.filter(report=self.discordance_report).select_related(
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
                if self.discordance_report.is_latest:
                    # classification is still outstanding, check to see if there pending changes
                    if flag := drc.classification_original.classification.flag_collection.get_open_flag_of_type(classification_flag_types.classification_pending_changes):
                        clinical_significance_to = flag.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY) or "unknown"
                        pending = True

                if not clinical_significance_to:
                    clinical_significance_to = drc.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

            group_counts[ClassificationLabSummaryEntry(
                lab=drc.classification_original.classification.lab,
                clinical_significance_from=clinical_significance_from,
                clinical_significance_to=clinical_significance_to,
                pending=pending
            )] += 1

        return sorted([ClassificationLabSummary(
            group=group,
            is_internal=group.lab in self.perspective.labs_if_not_admin,
            count=count
        ) for group, count in group_counts.items()])


@dataclass(frozen=True)
class DiscordanceReportSummaryCount:
    lab: Optional[Lab]
    count: int

    @property
    def is_internal(self):
        return self.lab is None

    def __lt__(self, other):
        # if we want internal to go first no matter the count
        if self.is_internal != other.is_internal:
            return self.is_internal
        if self.count != other.count:
            return self.count > other.count

        return self.lab < other.lab


@dataclass(frozen=True)
class DiscordanceReportTableData:
    perspective: LabPickerData
    summaries: List[DiscordanceReportRowData]
    inactive_summaries: List[DiscordanceReportRowData]

    counts: List[DiscordanceReportSummaryCount]

    def __len__(self):
        return len(self.summaries)

    def __bool__(self):
        return bool(self.summaries)

    def labs(self) -> List[Lab]:
        return sorted(self.perspective.your_labs)

    def labs_quick_str(self) -> str:
        if len(self.perspective.your_labs) == 1:
            lab = list(self.perspective.your_labs)[0]
            return str(lab)
        else:
            return "your assigned labs"

    @property
    def genome_build(self) -> GenomeBuild:
        return self.perspective.genome_build

    @staticmethod
    def create(perspective: LabPickerData, discordance_reports: Iterable[DiscordanceReport]) -> 'DiscordanceReportTableData':
        internal_count = 0
        by_lab: Dict[Lab, int] = defaultdict(int)
        summaries: List[DiscordanceReportRowData] = list()
        inactive_summaries: List[DiscordanceReportRowData] = list()

        for dr in discordance_reports:
            summary = DiscordanceReportRowData(discordance_report=dr, perspective=perspective)
            if summary.is_valid_including_withdraws:
                if summary.is_requiring_attention:
                    summaries.append(summary)
                    if summary.is_internal:
                        internal_count += 1
                    else:
                        for lab in summary.other_labs:
                            by_lab[lab] += 1
                else:
                    inactive_summaries.append(summary)

        counts = sorted([DiscordanceReportSummaryCount(lab=key, count=value) for key, value in by_lab.items()])
        if internal_count:
            counts.insert(0, DiscordanceReportSummaryCount(lab=None, count=internal_count))

        return DiscordanceReportTableData(perspective=perspective, summaries=summaries, inactive_summaries=inactive_summaries, counts=counts)


class DiscordanceAction:

    def __init__(self, code: int, short_description: str, long_description: str = None, no_longer_considered: bool = False):
        self.code = code
        self.short_description = short_description
        self.long_description = long_description
        self.no_longer_considered = no_longer_considered
        if not long_description:
            self.long_description = self.short_description

    def __eq__(self, other):
        if isinstance(other, DiscordanceAction):
            return self.code == other.code
        return False

    def __hash__(self):
        return self.code

    def __str__(self):
        return f'({self.code}) {self.short_description}'


DatedDiscordanceAction = typing.NamedTuple('DatedDiscordanceAction', [('action', DiscordanceAction), ('date', Optional[datetime])])


class DiscordanceActionsLog:

    NO_CHANGE = DiscordanceAction(0, "No change in classification")
    NOT_REVIEWED = DiscordanceAction(1, "Not reviewed")

    CHANGE_NO_REASON_GIVEN = DiscordanceAction(10, "Reclassified no reason given")
    CHANGE_UNKNOWN = DiscordanceAction(11, "Reclassified with reason (not supported in current code)")

    CHANGE_DD = DiscordanceAction(20, "Reclassified after discordance discussion")
    CHANGE_ND = DiscordanceAction(21, "Reclassified after summation of data")
    CHANGE_IR = DiscordanceAction(22, "Reclassified after internal review")

    CLINICAL_GROUP_CHANGED = DiscordanceAction(100, "Clinical grouping changed", no_longer_considered=True)
    CLASSIFICATION_WITHDRAWN = DiscordanceAction(101, "Classification withdrawn", no_longer_considered=True)

    ALL_ACTIONS = [
        NO_CHANGE,
        NOT_REVIEWED,
        CHANGE_NO_REASON_GIVEN,
        CHANGE_UNKNOWN,
        CHANGE_DD,
        CHANGE_ND,
        CHANGE_IR,
        CLINICAL_GROUP_CHANGED,
        CLASSIFICATION_WITHDRAWN
    ]

    def __init__(self):
        self.actions = list()
        self.internal_reviewed = None

    def add(self, action: DiscordanceAction, date: Optional[datetime] = None):
        self.actions.append(DatedDiscordanceAction(action=action, date=date))
        self.internal_reviewed = None

    @property
    def no_longer_considered(self):
        for dated_action in self.actions:
            if dated_action.action.no_longer_considered:
                return True
        return False

    def __str__(self):
        return f'internal review = {self.internal_reviewed}, actions = {self.actions}'


class DiscordanceReportClassificationRelationManager(models.Manager):

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.select_related('classification_original',
                                 'classification_final',
                                 'classification_final__classification',
                                 'classification_final__classification__lab',
                                 'classification_final__classification__lab__organization')


class DiscordanceReportClassification(TimeStampedModel):
    objects = DiscordanceReportClassificationRelationManager()
    report = models.ForeignKey(DiscordanceReport, on_delete=CASCADE)
    classification_original = models.ForeignKey(ClassificationModification, related_name='+', on_delete=CASCADE)
    classification_final = models.ForeignKey(ClassificationModification, related_name='+', on_delete=CASCADE, null=True, blank=True)
    clinical_context_final = models.ForeignKey(ClinicalContext, on_delete=PROTECT, null=True, blank=True)
    withdrawn_final = models.BooleanField(default=False)

    @lazy
    def classification_effective(self) -> ClassificationModification:
        """
        The final state of the classification (or the current if the report is still open)
        """
        if self.classification_final:
            return self.classification_final
        return ClassificationModification.objects.filter(
            classification=self.classification_original.classification,
            is_last_published=True
        ).select_related('classification', 'classification__lab').get()

    @lazy
    def clinical_context_effective(self) -> ClinicalContext:
        """
        The final state clinical context of the classification (or the current if the report is still open)
        """
        if self.clinical_context_final:
            return self.clinical_context_final
        return self.classification_original.classification.clinical_context

    @lazy
    def withdrawn_effective(self) -> bool:
        """
        The final "withdrawn" of the classification (or the current if the report is still open)
        """
        if self.classification_final:
            return self.withdrawn_final
        return self.classification_original.classification.withdrawn

    def close(self):
        self.clinical_context_final = self.clinical_context_effective
        self.withdrawn_final = self.withdrawn_effective
        self.classification_final = self.classification_effective
        self.save()

    @lazy
    def action_log(self) -> DiscordanceActionsLog:
        actions = DiscordanceActionsLog()

        if self.withdrawn_effective:
            # if withdrawn, nothing else matters
            actions.add(DiscordanceActionsLog.CLASSIFICATION_WITHDRAWN)
            return actions

        # Need to revisit this, how long do we give for flags to be updated after the discordance report has been closed
        # Currently 1 day, should it be 1 month? (but only if the relevant flag was opened during the discordance period
        # e.g. a significance_change flag opened on day 1, and finally filled in on day 20)
        flag_collection = self.classification_original.classification.flag_collection_safe
        start = self.report.created
        end = self.report.report_completed_date

        # find all flag comments from 1 minute before discordance was recorded
        # (as we want to get specific values)
        # but only for flags that are now closed.
        start -= timedelta(minutes=1)
        if end:
            # add configurable time to get stuff that happened after discordance was completed
            # gives time for people to fill in things
            end += timedelta(days=settings.DISCORDANCE_REPORT_LEEWAY)

        relevant_comments_qs = FlagComment.objects.filter(
            flag__collection=flag_collection,
            # resolution__status=FlagStatus.CLOSED,
            # flag__resolution__status=FlagStatus.CLOSED,
            flag__flag_type__in=[
                classification_flag_types.internal_review,
                classification_flag_types.significance_change,
            ],
        ).order_by('created')

        relevant_comments_qs = relevant_comments_qs.filter(created__gte=start)
        if end:
            relevant_comments_qs = relevant_comments_qs.filter(created__lte=end)

        # note if we had a significant change we treat it as if the record was reviewed
        # e.g. there's no
        had_significant_change = ClinicalSignificance.is_significant_change(
            old_classification=self.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
            new_classification=self.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        )

        if self.report.clinical_context != self.clinical_context_effective:
            actions.add(DiscordanceActionsLog.CLINICAL_GROUP_CHANGED)

        # TODO
        # flags = relevant_comments_qs.values_list('flag', flat=True).distinct() # go through in order and overwrite, so latest is most important

        has_reclass_reason = False
        internal_reviewed = False

        processed_flag = set()
        for flag_comment in relevant_comments_qs:  # : :type flag_comment: FlagComment
            # want to update resolution to be
            flag = flag_comment.flag  # : :type flag: Flag
            if flag.id in processed_flag:
                continue
            processed_flag.add(flag.id)

            if flag.flag_type == classification_flag_types.internal_review and flag.resolution.status == FlagStatus.CLOSED:
                internal_reviewed = True  # maybe we should get more info

            elif flag.flag_type == classification_flag_types.significance_change:
                resolution = flag.resolution
                has_reclass_reason = True
                if resolution and resolution.status == FlagStatus.CLOSED:
                    if resolution.id == 'sc_discordance_discussion':
                        actions.add(DiscordanceActionsLog.CHANGE_DD, flag_comment.created)
                    elif resolution.id == 'sc_new_data':
                        actions.add(DiscordanceActionsLog.CHANGE_ND, flag_comment.created)
                    elif resolution.id == 'sc_internal_review':
                        actions.add(DiscordanceActionsLog.CHANGE_IR, flag_comment.created)
                    else:
                        actions.add(DiscordanceActionsLog.CHANGE_UNKNOWN, flag_comment.created)
                else:
                    actions.add(DiscordanceActionsLog.CHANGE_NO_REASON_GIVEN, flag_comment.created)

        if internal_reviewed:
            actions.internal_reviewed = True

        if internal_reviewed and not had_significant_change:
            actions.add(DiscordanceActionsLog.NO_CHANGE)

        if not internal_reviewed and not had_significant_change:
            actions.add(DiscordanceActionsLog.NOT_REVIEWED)

        if had_significant_change and not has_reclass_reason:
            actions.add(DiscordanceActionsLog.CHANGE_NO_REASON_GIVEN)

        return actions
