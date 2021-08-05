from typing import Optional, List

from django.db import models
from django.db.models import QuerySet
from django.urls import reverse
from lazy import lazy
from model_utils.models import TimeStampedModel

from classification.json_utils import ValidatedJson, JsonObjType
from classification.models import ClassificationModification, ConditionResolved
from snpdb.models import ClinVarKey, Allele
from django.utils.timezone import now
import copy


class ClinVarAllele(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar allele"

    allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=models.CASCADE)

    classifications_missing_condition = models.IntegerField(default=0)
    submissions_valid = models.IntegerField(default=0)
    submissions_invalid = models.IntegerField(default=0)
    last_evaluated = models.DateTimeField(default=now)

    def __str__(self):
        return f"{self.allele} {self.clinvar_key}"


class ClinVarExportStatus(models.TextChoices):
    NEW_SUBMISSION = "N", "New Submission"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    CHANGES_PENDING = "C", "Changes Pending"
    UP_TO_DATE = "D", "Up to Date"
    IN_ERROR = "E", "Error"


class ClinVarReleaseStatus(models.TextChoices):
    WHEN_READY = "R", "Release When Ready"
    ON_HOLD = "H", "On Hold"


class ClinVarExport(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export"

    clinvar_allele = models.ForeignKey(ClinVarAllele, null=True, blank=True, on_delete=models.CASCADE)
    condition = models.JSONField()
    classification_based_on = models.ForeignKey(ClassificationModification, null=True, blank=True, on_delete=models.CASCADE)
    scv = models.TextField(null=True, blank=True)  # if not set yet
    status = models.CharField(max_length=1, choices=ClinVarExportStatus.choices, default=ClinVarExportStatus.NEW_SUBMISSION)
    release_status = models.CharField(max_length=1, choices=ClinVarReleaseStatus.choices, default=ClinVarReleaseStatus.WHEN_READY)
    last_evaluated = models.DateTimeField(default=now)
    submission_body_validated = models.JSONField(null=False, blank=False, default=dict)

    def get_absolute_url(self):
        return reverse('clinvar_export', kwargs={'pk': self.pk})

    def friendly_submission_status_str(self):
        if self.status == ClinVarExportStatus.NEW_SUBMISSION:
            return "New Submission"
        elif self.status == ClinVarExportStatus.CHANGES_PENDING:
            return "Valid (with changes pending)"
        elif self.status == ClinVarExportStatus.UP_TO_DATE:
            return "Valid (and up to date)"
        else:
            return "Submission has Errors"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @lazy
    def _condition_resolved(self) -> ConditionResolved:
        return ConditionResolved.from_dict(self.condition)

    @property
    def condition_resolved(self) -> ConditionResolved:
        return ConditionResolved.from_dict(self.condition)

    @condition_resolved.setter
    def condition_resolved(self, new_condition: ConditionResolved):
        self.condition = new_condition.to_json()
        lazy.invalidate(self, '_condition_resolved')

    def update_classification(self, new_classification_based_on: Optional[ClassificationModification]):
        if self.classification_based_on != new_classification_based_on:
            lazy.invalidate(self, 'submission_body')
            self.classification_based_on = new_classification_based_on
            self.update()

    @staticmethod
    def new_condition(clinvar_allele: ClinVarAllele, condition: ConditionResolved, candidate: Optional[ClassificationModification]) -> 'ClinVarExport':
        cc = ClinVarExport(clinvar_allele=clinvar_allele, condition=condition.to_json())
        cc.update_classification(candidate)
        return cc

    @lazy
    def submission_body(self) -> ValidatedJson:
        return ValidatedJson.deserialize(self.submission_body_validated)

    @property
    def submission_full(self) -> ValidatedJson:
        """
        Returns the data that should be directly copied into a ClinVarBatch
        """
        content = copy.deepcopy(self.submission_body)
        # ValidatedJson is meant to be immutable, but will make an exception here
        if scv := self.scv:
            content["recordStatus"] = "update"
            content["clinvarAccession"] = scv
        else:
            content["recordStatus"] = "novel"
        return content

    @property
    def last_submission(self) -> Optional['ClinVarExportSubmissionBatch']:
        return self.clinvarexportsubmission_set.exclude(submission_batch__status=ClinVarExportSubmissionBatchStatus.REJECTED).order_by('-created').first()

    def submission_body_previous(self) -> Optional[JsonObjType]:
        # ignore rejected submissions
        if last_submission := self.last_submission:
            return last_submission.submission_body
        return None

    def update(self):
        """
        Uses the linked data to newly generate JSON
        Will not update the classification used, update_classification should generally be called externally
        """
        from classification.models.clinvar_export_convertor import ClinVarExportConverter
        current_validated_json_body: ValidatedJson = ClinVarExportConverter(clinvar_export_record=self).as_json
        self.submission_body_validated = current_validated_json_body.serialize()
        lazy.invalidate(self, 'submission_body')

        status: ClinVarExportStatus
        if current_validated_json_body.has_errors:
            status = ClinVarExportStatus.IN_ERROR
        else:
            if previous_submission := self.submission_body_previous():
                if previous_submission != current_validated_json_body.pure_json():
                    status = ClinVarExportStatus.CHANGES_PENDING
                else:
                    status = ClinVarExportStatus.UP_TO_DATE
            else:
                status = ClinVarExportStatus.NEW_SUBMISSION
        self.status = status
        self.save()


class ClinVarExportSubmissionBatchStatus(models.TextChoices):
    AWAITING_UPLOAD = "W", "Awaiting Upload"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    UPLOADING = "U", "Uploading"
    SUBMITTED = "S", "Submitted"
    REJECTED = "R", "Rejected"


class ClinVarExportSubmissionBatch(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar Export submission batch"

    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.CASCADE)
    submission_version = models.IntegerField()
    status = models.CharField(max_length=1, choices=ClinVarExportSubmissionBatchStatus.choices, default=ClinVarExportSubmissionBatchStatus.AWAITING_UPLOAD)
    submitted_json = models.JSONField(null=True, blank=True)  # leave this as blank until we've actually uploaded data
    response_json = models.JSONField(null=True, blank=True)  # what did ClinVar return

    def get_absolute_url(self):
        return reverse('clinvar_export_batch', kwargs={'pk': self.pk})

    def __str__(self):
        return f"ClinVar Submission Batch : {self.id} - {self.get_status_display()}"

    def to_json(self) -> JsonObjType:
        return {
            "behalfOrgId": self.clinvar_key.behalf_org_id,
            "clinVarSubmissions": [submission.submission_full for submission in self.clinvarexportsubmission_set.order_by('created')],
            "submissionName": f"submission_{self.id}"
        }

    @staticmethod
    def create_batches(qs: QuerySet) -> List['ClinVarExportSubmissionBatch']:
        all_batches: List[ClinVarExportSubmissionBatch] = list()

        current_batch: Optional[ClinVarExportSubmissionBatch] = None
        current_batch_size = 0

        qs = qs.order_by('clinvar_allele__clinvar_key')
        qs = qs.filter(status__in=[ClinVarExportStatus.NEW_SUBMISSION, ClinVarExportStatus.CHANGES_PENDING])
        qs = qs.select_related('clinvar_allele', 'clinvar_allele__clinvar_key')
        record: ClinVarExport
        for record in qs:
            full_current = record.submission_full
            if not full_current.has_errors:  # should never have errors as we're filtering on new submission changes pending

                if current_batch_size == 0 or current_batch_size == 10000 or current_batch.clinvar_key != record.clinvar_allele.clinvar_key:
                    # if current_batch is not None:
                    #    current_batch.refresh_from_db()

                    current_batch = ClinVarExportSubmissionBatch(clinvar_key=record.clinvar_allele.clinvar_key,
                                                                 submission_version=1)
                    current_batch.save()
                    current_batch_size = 0
                    all_batches.append(current_batch)

                ClinVarExportSubmission(
                    clinvar_export=record,
                    classification_based_on=record.classification_based_on,
                    submission_batch=current_batch,
                    submission_full=full_current.pure_json(),
                    submission_body=record.submission_body.pure_json(),
                    submission_version=1
                ).save()
                record.status = ClinVarExportStatus.UP_TO_DATE
                record.save()

        return all_batches


class ClinVarExportSubmission(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export submission"

    clinvar_export = models.ForeignKey(ClinVarExport, on_delete=models.CASCADE)  # if there's been an actual submission, don't allow deletes
    classification_based_on = models.ForeignKey(ClassificationModification, on_delete=models.PROTECT)
    submission_batch = models.ForeignKey(ClinVarExportSubmissionBatch, on_delete=models.CASCADE)
    submission_full = models.JSONField()  # the full data included in the batch submission

    submission_body = models.JSONField()  # used to see if there are any changes since last submission (other than going from novel to update)
    submission_version = models.IntegerField()


"""
@receiver(post_save, sender=ClinVarExport)
def set_condition_alias_permissions(sender, created: bool, instance: ClinVarExport, **kwargs):  # pylint: disable=unused-argument
    if created:
        group = instance.lab.group
        assign_perm(ClinVarExport.get_read_perm(), group, instance)
        assign_perm(ClinVarExport.get_write_perm(), group, instance)


class ClinVarExportSubmission(TimeStampedModel, GuardianPermissionsMixin):
    clinvar_export = models.ForeignKey(ClinVarExport, on_delete=CASCADE)
    evidence = models.JSONField()
    submission_status = models.TextField()


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):

    cve: ClinVarExport
    if cve := ClinVarExport.objects.filter(classification_based_on__classification=classification).first():
        cve.update_with(newly_published)
        cve.save()


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):  # pylint: disable=unused-argument
    # Keeps condition_text_match in sync with the classifications when withdraws happen/finish
    flag = flag_comment.flag
    if flag.flag_type == flag_types.classification_flag_types.classification_withdrawn:
        cl: Classification
        if cl := Classification.objects.filter(flag_collection=flag.collection.id).first():
            cve: ClinVarExport
            if cve := ClinVarExport.objects.filter(classification_based_on__classification=cl).first():
                cve.withdrawn = flag_comment.resolution.status == FlagStatus.OPEN
                cve.save()
"""