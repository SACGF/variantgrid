import html
from enum import Enum
from functools import cached_property
from typing import Optional, Any, Union

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User, Group
from django.core.exceptions import PermissionDenied
from django.db import models, transaction
from django.db.models import TextChoices
from django.db.models.deletion import PROTECT, CASCADE
from django.urls.base import reverse
from django.utils import timezone
from django.utils.safestring import mark_safe
from django_extensions.db.models import TimeStampedModel

from classification.enums.classification_enums import SpecialEKeys
from classification.enums.discordance_enums import DiscordanceReportResolution, ContinuedDiscordanceReason
from classification.models.classification import ClassificationModification, Classification
from classification.models.clinical_context_models import ClinicalContext, ClinicalContextRecalcTrigger
from classification.models.clinical_context_models import ClinicalContextChangeData
from classification.models.flag_types import classification_flag_types
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from library.preview_request import PreviewModelMixin, PreviewKeyValue, PreviewData
from library.utils import invalidate_cached_property
from library.utils.django_utils import refresh_for_update
from review.models import ReviewableModelMixin, Review
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab, GenomeBuild


discordance_change_signal = django.dispatch.Signal()  # args: "discordance_report", "clinical_context_change_data:ClinicalContextChangeData"


class NotifyLevel(str, Enum):
    NEVER_NOTIFY = "never-notify"
    NOTIFY_IF_CHANGE = "notify-if-change"
    ALWAYS_NOTIFY = "always-notify"


class DiscordanceReport(TimeStampedModel, ReviewableModelMixin, PreviewModelMixin):

    resolution = models.TextField(default=DiscordanceReportResolution.ONGOING, choices=DiscordanceReportResolution.CHOICES, max_length=1, null=True, blank=True)
    # TODO remove continued discordance reason, it should be redundant to notes
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

    admin_note = models.TextField(null=False, blank=True, default='')

    @classmethod
    def preview_category(cls) -> str:
        return "Discordance Report"

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-arrow-down-up-across-line"

    @classmethod
    def preview_enabled(cls) -> bool:
        return settings.DISCORDANCE_ENABLED

    @property
    def preview(self) -> 'PreviewData':

        from classification.views.discordance_report_views import DiscordanceReportTemplateData
        drtd = DiscordanceReportTemplateData(self.pk, user=admin_bot())

        c_hgvs_key_values = []
        for c_hgvs in drtd.c_hgvses:
            c_hgvs_key_values.append(
                PreviewKeyValue(key=f"{c_hgvs.genome_build} c.HGVS", value=str(c_hgvs), dedicated_row=True)
            )

        status_text: str
        if self.is_pending_concordance and self.is_latest:
            status_text = "Pending Concordance"
        else:
            status_text = self.get_resolution_display() or 'Active Discordance'

        # note there's also preview_extra_signal that provides the lab data
        return self.preview_with(
            identifier=f"DR_{self.pk}",
            summary_extra=
                [PreviewKeyValue(key="Status", value=status_text, dedicated_row=True)] +
                [PreviewKeyValue(key="Allele", value=f"{self.clinical_context.allele:CA}", dedicated_row=True)] +
                c_hgvs_key_values
        )

    class LabInvolvement(int, Enum):
        WITHDRAWN = 1
        ACTIVE = 2

    def get_absolute_url(self):
        return reverse('discordance_report', kwargs={"discordance_report_id": self.pk})

    def __str__(self):
        return f"({self.pk}) Discordance Report ({self.clinical_context.allele:CA}) {self.get_resolution_display() or 'Discordant'}"

    @property
    def metrics_logging_key(self) -> tuple[str, Any]:
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
    def close(self, clinical_context_change_data: ClinicalContextChangeData, expected_resolution: Optional[str] = None):
        """
        @param expected_resolution will be calculated, but if provided a ValueError will be raised if it doesn't match
        @param clinical_context_change_data reason this discordance is being closed
        """
        if self.clinical_context.discordance_status.is_discordant:
            self.resolution = DiscordanceReportResolution.CONTINUED_DISCORDANCE
        else:
            self.resolution = DiscordanceReportResolution.CONCORDANT

        if expected_resolution and expected_resolution != self.resolution:
            raise ValueError(f'Expected to close Discordance Report as {expected_resolution} but was {self.resolution}')

        self.resolved_text = clinical_context_change_data.cause_text
        self.report_completed_date = timezone.now()
        self.save()

        for drc in DiscordanceReportClassification.objects.filter(report=self):  # type: DiscordanceReportClassification
            drc.close()

        discordance_change_signal.send(DiscordanceReport, discordance_report=self, clinical_context_change_data=clinical_context_change_data)

    @transaction.atomic
    def update(self, clinical_context_change_data: ClinicalContextChangeData, notify_level: NotifyLevel = NotifyLevel.NOTIFY_IF_CHANGE):
        if not self.is_active:
            raise ValueError('Cannot update non active Discordance Report')

        if not self.id:
            self.save()

        existing_vms = set()
        existing_labs = set()
        for drc in DiscordanceReportClassification.objects.filter(report=self):
            existing_vms.add(drc.classification_original.classification_id)
            existing_labs.add(drc.classification_original.classification.lab)

        newly_added_labs: set[Lab] = set()
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

        invalidate_cached_property(self, 'discordance_report_classifications')
        invalidate_cached_property(self, 'involved_labs')

        # contents of existing_vms are now presumably records that changed clinical groupings

        if not self.clinical_context.discordance_status.is_discordant:
            # hey we've reached concordance, let's close this
            # close will fire off a notification event
            self.close(clinical_context_change_data=clinical_context_change_data, expected_resolution=DiscordanceReportResolution.CONCORDANT)
        else:
            if notify_level == NotifyLevel.NEVER_NOTIFY:
                pass
            elif notify_level == NotifyLevel.ALWAYS_NOTIFY:
                discordance_change_signal.send(DiscordanceReport, discordance_report=self, clinical_context_change_data=clinical_context_change_data)
            elif newly_added_labs:  # change is significant
                newly_added_labs_str = ", ".join(str(lab) for lab in newly_added_labs)
                clinical_context_change_data_cause = ClinicalContextChangeData(cause_text=f"{clinical_context_change_data.cause_text} and newly added labs {newly_added_labs_str}", cause_code=clinical_context_change_data.cause_code)
                discordance_change_signal.send(DiscordanceReport, discordance_report=self, clinical_context_change_data=clinical_context_change_data_cause)
            else:
                # the complete withdrawl of 1 lab means we might still want to close off triages
                discordance_change_signal.send(
                    DiscordanceReport,
                    discordance_report=self,
                    clinical_context_change_data=clinical_context_change_data.with_notify_worthy(notify=False)
                )

    @property
    def all_actively_involved_labs(self) -> set[Lab]:
        return {lab for lab, status in self.involved_labs.items() if status == DiscordanceReport.LabInvolvement.ACTIVE}

    @property
    def reviewing_labs(self) -> set[Lab]:
        return self.all_actively_involved_labs

    def post_review_url(self, review: Review) -> str:
        return reverse('discordance_report_review_action', kwargs={"review_id": review.pk})

    @property
    def is_review_locked(self) -> bool:
        if self.resolution:
            return True

        if self != DiscordanceReport.objects.filter(clinical_context=self.clinical_context).order_by('-report_started_date').first():
            return True

        return False

    @cached_property
    def discordance_report_classifications(self) -> list['DiscordanceReportClassification']:
        return list(self.discordancereportclassification_set.select_related(
            'classification_original__classification__clinical_context',
            'classification_final__classification',
            'classification_original__classification__allele_info__grch37',
            'classification_original__classification__allele_info__grch38',
            'classification_final__classification__allele_info__grch37',
            'classification_final__classification__allele_info__grch38'
        ).all())

    @cached_property
    def involved_labs(self) -> dict[Lab, LabInvolvement]:
        lab_status: dict[Lab, DiscordanceReport.LabInvolvement] = {}
        for drc in self.discordance_report_classifications:
            effective_c: Classification = drc.classification_effective.classification
            lab = effective_c.lab
            status = DiscordanceReport.LabInvolvement.WITHDRAWN if effective_c.withdrawn else DiscordanceReport.LabInvolvement.ACTIVE
            lab_status[lab] = max(lab_status.get(lab, 0), status)
        return lab_status

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        return self.user_is_involved(user_or_group)

    def check_can_view(self, user_or_group: Union[User, Group]):
        if not self.can_view(user_or_group):
            msg = f"You do not have READ permission to view {self.pk}"
            raise PermissionDenied(msg)

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
        ccd = ClinicalContextChangeData(
            cause_text=cause,
            cause_code=ClinicalContextRecalcTrigger.RE_OPEN
        )

        report.update(clinical_context_change_data=ccd, notify_level=NotifyLevel.ALWAYS_NOTIFY)
        return report

    @property
    def has_significance_changed(self):
        existing_vms = {}
        for dr in self.discordance_report_classifications:
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
            # what was once "continued discordance" is concordant, re-open, so we can instantly close it
            return True

        existing_labs: set[Lab] = set()
        for drc in self.discordance_report_classifications:
            existing_labs.add(drc.classification_original.classification.lab)

        # a new lab has submitted since continued discordance
        for cc in self.clinical_context.classifications_qs:
            if cc.lab not in existing_labs:
                return True
        return False

    @staticmethod
    def latest_report(clinical_context: ClinicalContext) -> 'DiscordanceReport':
        return clinical_context.discordancereport_set.order_by('-created').first()

    @property
    def is_latest(self):
        return DiscordanceReport.latest_report(self.clinical_context) == self

    @staticmethod
    def update_latest(clinical_context: ClinicalContext, clinical_context_change_data: ClinicalContextChangeData, update_flags: bool):
        try:
            latest_report = DiscordanceReport.latest_report(clinical_context=clinical_context)
            if latest_report:
                if latest_report.is_active:
                    latest_report.update(clinical_context_change_data=clinical_context_change_data)
                    return latest_report
                if latest_report.resolution == DiscordanceReportResolution.CONCORDANT:
                    # most recent report ended with discordance, can start fresh
                    pass
                elif latest_report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
                    # create a new report if data has changed significantly
                    if latest_report.should_reopen_continued_discordance:
                        return latest_report.reopen_continued_discordance(cause=f'Change after previous report marked as continued discordance - {clinical_context_change_data.cause_text}')
                    return None

            if clinical_context.is_discordant():
                report = DiscordanceReport(clinical_context=clinical_context, cause_text=clinical_context_change_data.cause_text)
                report.update(clinical_context_change_data=clinical_context_change_data)
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
        discordant_classifications: set[int] = set()
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

    def actively_discordant_classification_ids(self) -> set[int]:
        if not self.is_important:
            return set()
        classifications = set()
        for dr in self.discordance_report_classifications:
            if dr.clinical_context_effective == self.clinical_context and not dr.withdrawn_effective:
                classifications.add(dr.classification_original.classification.id)
        return classifications

    @property
    def all_classification_modifications(self) -> list[ClassificationModification]:
        return [drc.classification_effective for drc in self.discordance_report_classifications]

    @cached_property
    def is_medically_significant(self):
        ds = self.clinical_context.discordance_status

        if ds.is_discordant and not ds.pending_concordance:
            for cm in self.all_classification_modifications:
                if not cm.classification.withdrawn:
                    if "P" in cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                        return True
        return False

    def all_c_hgvs(self, genome_build: Optional[GenomeBuild] = None) -> list[CHGVS]:
        if not genome_build:
            genome_build = GenomeBuildManager.get_current_genome_build()
        c_hgvs = set()
        for cm in self.all_classification_modifications:
            c_hgvs.add(cm.c_hgvs_best(genome_build))
        return sorted(c_hgvs)

    @property
    def is_pending_concordance(self):
        return self.clinical_context.discordance_status.pending_concordance


# TODO all the below classes are utilities, consider moving them out

class DiscordanceReportNextStep(str, Enum):
    UNANIMOUSLy_COMPLEX = "C"
    AWAITING_YOUR_TRIAGE = "T"
    AWAITING_YOUR_AMEND = "A"
    AWAITING_OTHER_LAB = "O"
    TO_DISCUSS = "D"


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

    @cached_property
    def classification_effective(self) -> ClassificationModification:
        """
        The final state of the classification (or the current if the report is still open)
        """
        if self.classification_final:
            return self.classification_final
        return ClassificationModification.objects.filter(
            classification=self.classification_original.classification,
            is_last_published=True
        ).select_related('classification', 'classification__lab', 'classification__lab__organization').get()

    @cached_property
    def clinical_context_effective(self) -> ClinicalContext:
        """
        The final state clinical context of the classification (or the current if the report is still open)
        """
        if self.clinical_context_final:
            return self.clinical_context_final
        return self.classification_original.classification.clinical_context

    @cached_property
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


class DiscordanceReportTriageStatus(TextChoices):
    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Classification"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"


class DiscordanceReportTriage(PreviewModelMixin, TimeStampedModel):
    discordance_report = models.ForeignKey(DiscordanceReport, on_delete=CASCADE)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    triage_status = models.TextField(max_length=1, choices=DiscordanceReportTriageStatus.choices, default=DiscordanceReportTriageStatus.PENDING)
    note = models.TextField(null=True, blank=True)
    triage_date = models.DateField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, blank=True, on_delete=PROTECT)
    closed = models.BooleanField(default=False)
    """
    These tags are more dynamic, so keep it configurable in JSON rather than a database schema for now
    """

    def can_write(self, user: User) -> bool:
        return self.lab.is_member(user, admin_check=True) and not self.closed

    @property
    def is_outstanding(self):
        # another example of pending_concordance neeing to be an official status
        return self.triage_status == DiscordanceReportTriageStatus.PENDING and \
            self.discordance_report.clinical_context.discordance_status.is_discordant and \
            not self.discordance_report.is_pending_concordance

    @staticmethod
    def outstanding_for_lab(lab: Lab):
        return (drt for drt in DiscordanceReportTriage.objects.filter(lab=lab, review_status=DiscordanceReportTriageStatus.PENDING) if drt.is_outstanding)

    @cached_property
    def note_html(self):
        parts = []
        if user := self.user:
            parts.append(f"Updated by: { user }")
        if triage_date := self.triage_date:
            parts.append(f"On: { triage_date }")
        parts.append(f"To: {self.get_triage_status_display()}")
        if note := self.note:
            parts.append("Note: " + html.escape(note, quote=True))
        if parts:
            return "<br/>".join(parts)

    @classmethod
    def preview_category(cls) -> str:
        return "Discordance Report Triage"

    def get_absolute_url(self):
        return reverse('discordance_report', kwargs={"discordance_report_id": self.discordance_report.pk})

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-arrow-down-up-across-line"

    @property
    def preview(self) -> 'PreviewData':
        title = f"DR_ {self.discordance_report.pk}"

        extras = [
            PreviewKeyValue(key="Status", value=self.get_triage_status_display()),
            PreviewKeyValue(key="Modified", value=self.modified),
            PreviewKeyValue(key="Lab", value=self.lab),
        ]

        return PreviewData.for_object(
            obj=self,
            category="Discordance Report Triage",
            title=title,
            icon=self.preview_icon(),
            internal_url=self.get_absolute_url(),
            summary_extra=extras
        )

    class Meta:
        unique_together = ('discordance_report', 'lab')


def ensure_discordance_report_triages_for(dr: DiscordanceReport):
    with transaction.atomic():
        dr = refresh_for_update(dr)

        resolved = dr.resolution is not None or dr.is_pending_concordance

        if not resolved:
            all_actively_involved_labs = dr.all_actively_involved_labs

            for lab in all_actively_involved_labs:
                DiscordanceReportTriage.objects.get_or_create({}, discordance_report=dr, lab=lab)

            DiscordanceReportTriage.objects.filter(discordance_report=dr).filter(lab__in=all_actively_involved_labs).update(closed=False)
            # the below will happen if a lab has withdrawn out of a discordance
            DiscordanceReportTriage.objects.filter(discordance_report=dr).exclude(lab__in=all_actively_involved_labs).update(closed=True)
        else:
            # it is resolved, everything should be closed
            DiscordanceReportTriage.objects.filter(discordance_report=dr).update(closed=True)


def ensure_discordance_report_triages_bulk():
    for dr in DiscordanceReport.objects.all():
        ensure_discordance_report_triages_for(dr)
