from typing import Optional, List, Iterable

from django.db import models, transaction
from django.db.models import QuerySet, TextChoices
from django.urls import reverse
from lazy import lazy
from model_utils.models import TimeStampedModel

from uicore.json.validated_json import ValidatedJson
from uicore.json.json_types import JsonObjType
from classification.models import ClassificationModification, ConditionResolved
from snpdb.models import ClinVarKey, Allele
from django.utils.timezone import now
import copy


CLINVAR_EXPORT_CONVERSION_VERSION = 2


class ClinVarAllele(TimeStampedModel):
    """
    Wraps an Allele with a ClinVarKey to make it easier to keep track of our submissions
    (Provides a little bit of bloat, but does give us allele level data)
    """

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
    """
    The status of an Export as compared to the last batch.
    Note doesn't involve itself in the status as data progresses
    """

    NEW_SUBMISSION = "N", "New Submission"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    CHANGES_PENDING = "C", "Changes Pending"
    UP_TO_DATE = "D", "Up to Date"
    IN_ERROR = "E", "Error"


class ClinVarReleaseStatus(models.TextChoices):
    """
    As determined by the user currently
    """
    WHEN_READY = "R", "Release When Ready"
    ON_HOLD = "H", "On Hold"


class ClinVarExport(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export"

    clinvar_allele = models.ForeignKey(ClinVarAllele, null=True, blank=True, on_delete=models.CASCADE)
    condition = models.JSONField()
    classification_based_on = models.ForeignKey(ClassificationModification, null=True, blank=True, on_delete=models.CASCADE)
    scv = models.TextField(blank=True, default='')
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
        cc = ClinVarExport(
            clinvar_allele=clinvar_allele,
            condition=condition.to_json()
        )
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
    def previous_submission(self) -> Optional['ClinVarExportSubmission']:
        return self.clinvarexportsubmission_set.exclude(submission_batch__status=ClinVarExportBatchStatus.REJECTED).order_by('-created').first()

    def update(self):
        """
        Uses the linked data to newly generate JSON
        Will not update the classification used, update_classification should generally be called externally
        """
        from classification.models.clinvar_export_convertor import ClinVarExportConverter

        current_validated_json_body: ValidatedJson = ClinVarExportConverter(clinvar_export_record=self).as_validated_json
        self.submission_body_validated = current_validated_json_body.serialize()
        lazy.invalidate(self, 'submission_body')

        status: ClinVarExportStatus
        if current_validated_json_body.has_errors:
            status = ClinVarExportStatus.IN_ERROR
        else:
            if previous_submission := self.previous_submission:
                if previous_submission.submission_body != current_validated_json_body.pure_json():
                    status = ClinVarExportStatus.CHANGES_PENDING
                else:
                    status = ClinVarExportStatus.UP_TO_DATE
            else:
                status = ClinVarExportStatus.NEW_SUBMISSION
        self.status = status
        self.save()


class ClinVarExportBatchStatus(models.TextChoices):
    """
    Overall status of a submission batch, talking to ClinVar takes a few steps
    """
    AWAITING_UPLOAD = "W", "Awaiting Upload"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    UPLOADING = "P", "Processing"
    SUBMITTED = "S", "Submitted"
    REJECTED = "R", "Rejected"


class ClinVarExportBatch(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export batch"

    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.CASCADE)
    submission_version = models.IntegerField()
    status = models.CharField(max_length=1, choices=ClinVarExportBatchStatus.choices, default=ClinVarExportBatchStatus.AWAITING_UPLOAD)
    submission_identifier = models.TextField(null=True, blank=True)
    file_url = models.TextField(null=True, blank=True)

    def __str__(self):
        return f"ClinVar Submission Batch : {self.id} - {self.get_status_display()}"

    def to_json(self) -> JsonObjType:
        return {
            "actions": [{
                "type": "AddData",
                "targetDb": "clinvar",
                "data": {
                    "content": {
                        "behalfOrgID": self.clinvar_key.behalf_org_id or "testorg",
                        "clinvarSubmission": [submission.submission_full for submission in
                                              self.clinvarexportsubmission_set.order_by('created')],
                        "submissionName": f"submission_{self.id}"
                    }
                }
            }]
        }

    def requests(self) -> Iterable['ClinVarExportRequest']:
        return self.clinvarexportrequest_set.order_by('created')

    @transaction.atomic()
    def reject(self):
        self.status = ClinVarExportBatchStatus.REJECTED
        self.save()
        exports = self.clinvarexportsubmission_set.values_list('clinvar_export', flat=True).distinct()
        for export in ClinVarExport.objects.filter(pk__in=exports):
            export.update()

    @staticmethod
    @transaction.atomic()
    def create_batches(qs: QuerySet[ClinVarExport], force_update: bool = False) -> List['ClinVarExportBatch']:
        all_batches: List[ClinVarExportBatch] = list()

        current_batch: Optional[ClinVarExportBatch] = None
        current_batch_size = 0

        qs = qs.exclude(release_status=ClinVarReleaseStatus.ON_HOLD)
        qs = qs.order_by('clinvar_allele__clinvar_key')

        if not force_update:
            qs = qs.filter(status__in=[ClinVarExportStatus.NEW_SUBMISSION, ClinVarExportStatus.CHANGES_PENDING])
        qs = qs.select_related('clinvar_allele', 'clinvar_allele__clinvar_key')
        record: ClinVarExport
        for record in qs:
            if force_update:
                record.update()
                # only have to do a check if we're not previously doing the filter
                if record.status in {ClinVarExportStatus.UP_TO_DATE, ClinVarExportStatus.IN_ERROR}:
                    continue

            full_current = record.submission_full
            if not full_current.has_errors:  # should never have errors as we're filtering on new submission changes pending

                if current_batch_size == 0 or current_batch_size == 10000 or current_batch.clinvar_key != record.clinvar_allele.clinvar_key:
                    # if current_batch is not None:
                    #    current_batch.refresh_from_db()

                    current_batch = ClinVarExportBatch(clinvar_key=record.clinvar_allele.clinvar_key,
                                                       submission_version=CLINVAR_EXPORT_CONVERSION_VERSION)
                    current_batch.save()
                    current_batch_size = 0
                    all_batches.append(current_batch)

                pure_json = record.submission_body.pure_json()
                local_id = pure_json["localID"]
                local_key = pure_json["localKey"]

                ClinVarExportSubmission(
                    clinvar_export=record,
                    classification_based_on=record.classification_based_on,
                    submission_batch=current_batch,
                    submission_full=full_current.pure_json(),
                    submission_body=record.submission_body.pure_json(),
                    submission_version=CLINVAR_EXPORT_CONVERSION_VERSION,
                    localId=local_id,
                    localKey=local_key
                ).save()
                current_batch_size += 1
                record.status = ClinVarExportStatus.UP_TO_DATE
                record.save()

        return all_batches


class ClinVarExportSubmissionStatus(TextChoices):
    WAITING = "W", "Waiting"
    SUCCESS = "S", "Success"
    ERROR = "E", "Error"


class ClinVarExportSubmission(TimeStampedModel):
    """
    Represents the content of a single ClinVarExport at a given point in time (within a submission batch)
    """

    class Meta:
        verbose_name = "ClinVar export submission"
        # below should be uncommented, but have to migrate any existing data first
        # unique_together = ("clinvar_export", "localKey")

    clinvar_export = models.ForeignKey(ClinVarExport, on_delete=models.CASCADE)  # if there's been an actual submission, don't allow deletes
    classification_based_on = models.ForeignKey(ClassificationModification, on_delete=models.PROTECT)
    submission_batch = models.ForeignKey(ClinVarExportBatch, on_delete=models.CASCADE)
    submission_full = models.JSONField()  # the full data included in the batch submission

    submission_body = models.JSONField()  # used to see if there are any changes since last submission (other than going from novel to update)
    submission_version = models.IntegerField()

    localId = models.TextField()  # should be static
    localKey = models.TextField()  # can change as the selected classification changes

    # individual record failure, batch can be in error too
    status = models.CharField(max_length=1, choices=ClinVarExportSubmissionStatus.choices, default=ClinVarExportSubmissionStatus.WAITING)
    scv = models.TextField(null=True, blank=True)
    response_json = models.JSONField(null=True, blank=True)


class ClinVarExportRequestType(TextChoices):
    INITIAL_SUBMISSION = "S", "Initial Submission"
    POLLING_SUBMISSION = "P", "Polling"
    RESPONSE_FILES = "F", "Response Files"


class ClinVarExportRequest(TimeStampedModel):
    """
    Represents a single request we sent to ClinVar
    """
    submission_batch = models.ForeignKey(ClinVarExportBatch, on_delete=models.CASCADE)
    request_type = models.CharField(max_length=1, choices=ClinVarExportRequestType.choices)
    url = models.TextField()
    request_json = models.JSONField(null=True, blank=True)  # some requests
    response_status_code = models.IntegerField()
    response_json = models.JSONField(null=True, blank=True)
    handled = models.BooleanField(default=False)


"""
# CODE for immediate syncing submissions to Exports

@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):
    pass


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):  # pylint: disable=unused-argument
    pass
"""