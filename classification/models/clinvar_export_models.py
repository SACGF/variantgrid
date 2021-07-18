import typing
from typing import Optional, List

from django.contrib.auth.models import User
from django.db import models
from django.db.models import QuerySet
from lazy import lazy
from model_utils.models import TimeStampedModel
from classification.json_utils import ValidatedJson, JsonObjType
from classification.models import ClassificationModification, ConditionResolved
from snpdb.models import ClinVarKey, Allele, Lab
import copy


class ClinVarAllele(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar allele"

    allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=models.CASCADE)

    def __str__(self):
        return f"{self.allele} {self.clinvar_key}"

    @staticmethod
    def clinvar_keys_for_user(user: User) -> QuerySet:
        """
        Ideally this would be on ClinVarKey but can't be due to ordering
        """
        if user.is_superuser:
            return ClinVarKey.objects.all()
        return ClinVarKey.objects.filter(pk__in=Lab.valid_labs_qs(user).filter(clinvar_key__isnull=False).select_related('clinvar_key'))


class ClinVarStatus(models.TextChoices):
    NEW_SUBMISSION = "N"  # new submission and changes pending often work the same, but might be useful to see at a glance, useful if we do approvals
    CHANGES_PENDING = "C"
    UP_TO_DATE = "D"
    IN_ERROR = "E"


class ClinVarExportRecord(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export record"

    clinvar_allele = models.ForeignKey(ClinVarAllele, null=True, blank=True, on_delete=models.CASCADE)
    condition = models.JSONField()
    classification_based_on = models.ForeignKey(ClassificationModification, null=True, blank=True, on_delete=models.CASCADE)
    scv = models.TextField(null=True, blank=True)  # if not set yet
    status = models.CharField(max_length=1, choices=ClinVarStatus.choices, default=ClinVarStatus.NEW_SUBMISSION)

    def __init__(self, *args, **kwargs):
        super(TimeStampedModel, self).__init__(*args, **kwargs)
        self.cached_condition: Optional[ConditionResolved] = None

    @property
    def condition_resolved(self) -> ConditionResolved:
        if not self.cached_condition:
            self.cached_condition = ConditionResolved.from_dict(self.condition)
        return self.cached_condition

    @condition_resolved.setter
    def condition_resolved(self, new_condition: ConditionResolved):
        self.condition = new_condition.as_json_minimal()
        self.cached_condition = new_condition

    def update_classification(self, new_classification_based_on: Optional[ClassificationModification]):
        if self.classification_based_on != new_classification_based_on:
            lazy.invalidate(self, 'submission_body_current')
            self.classification_based_on = new_classification_based_on
            self.status = self.calculate_status()
            # FIXME, check to see if we changed since last submission
            self.save()

    @staticmethod
    def new_condition(clinvar_allele: ClinVarAllele, condition: ConditionResolved, candidate: Optional[ClassificationModification]) -> 'ClinVarExportRecord':
        cc = ClinVarExportRecord(clinvar_allele=clinvar_allele, condition=condition.as_json_minimal())
        cc.update_classification(candidate)
        return cc

    @lazy
    def submission_body_current(self) -> ValidatedJson:
        """
        Returns the body of the submission, which wont change by going from novel -> update
        """
        from classification.models.clinvar_export_convertor import ClinVarExportConverter
        return ClinVarExportConverter(clinvar_export_record=self).as_json

    @property
    def submission_full_current(self) -> ValidatedJson:
        """
        Returns the data that should be directly copied into a ClinVarBatch
        """
        content = copy.deepcopy(self.submission_body_current)
        if scv := self.scv:
            content["recordStatus"] = "update"
            content["clinvarAccession"] = scv
        else:
            content["recordStatus"] = "novel"
        return content

    def submission_body_previous(self) -> Optional[JsonObjType]:
        if last_submission := self.clinvarexportrecordsubmission_set.order_by('-created').first():
            return last_submission.submission_body
        return None

    def calculate_status(self) -> ClinVarStatus:
        current_body = self.submission_body_current
        if current_body.has_errors:
            return ClinVarStatus.IN_ERROR
        else:
            if previous_submission := self.submission_body_previous():
                if previous_submission != self.submission_body_current.pure_json():
                    return ClinVarStatus.CHANGES_PENDING
                else:
                    return ClinVarStatus.UP_TO_DATE
            else:
                # no previous submission
                return ClinVarStatus.NEW_SUBMISSION


class ClinVarSubmissionBatchBuilder:

    def __init__(self, clinvar_key: ClinVarKey):
        self.clinvar_key = clinvar_key
        self.json = list()


class ClinVarSubmissionBatch(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar submission batch"

    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.PROTECT)
    submission_version = models.IntegerField()
    pass  # TODO add a bunch more fields when we know what they are

    def to_json(self) -> JsonObjType:
        return {
            "behalfOrgId": self.clinvar_key.behalf_org_id,
            "clinVarSubmissions": [submission.submission_full for submission in self.clinvarexportrecordsubmission_set.order_by('created')],
            "submissionName": f"submission_{self.id}"
        }

    @staticmethod
    def create_batches(self, qs: QuerySet) -> List['ClinVarSubmissionBatch']:
        all_batches: List[ClinVarSubmissionBatch] = list()

        current_batch: Optional[ClinVarSubmissionBatch] = None
        current_batch_size = 0

        qs = qs.order_by('clinvar_allele__clinvar_key')
        qs = qs.filter(status__in=[ClinVarStatus.NEW_SUBMISSION, ClinVarStatus.CHANGES_PENDING])
        qs = qs.select_related('clinvar_allele', 'clinvar_allele__clinvar_key')
        record: ClinVarExportRecord
        for record in qs:
            if current_batch_size == 0 or current_batch_size == 10000 or current_batch.clinvar_key != record.clinvar_allele.clinvar_key:
                current_batch = ClinVarSubmissionBatch(clinvar_key=record.clinvar_allele.clinvar_key, submission_version=1)
                current_batch_size = 0
                all_batches.append(current_batch)

            full_current = record.submission_full_current
            if not full_current.has_errors:
                ClinVarExportRecordSubmission.objects.create(
                    clinvar_candidate=record,
                    classification_based_on=record.classification_based_on,
                    submission_batch=current_batch,
                    submission_full=full_current.pure_json(),
                    submission_body=record.submission_body_current.pure_json(),
                    submission_version=1
                )
        return all_batches


class ClinVarExportRecordSubmission(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export record submission"

    clinvar_candidate = models.ForeignKey(ClinVarExportRecord, on_delete=models.PROTECT)  # if there's been an actual submission, don't allow deletes
    classification_based_on = models.ForeignKey(ClassificationModification, on_delete=models.PROTECT)
    submission_batch = models.ForeignKey(ClinVarSubmissionBatch, null=True, blank=True, on_delete=models.SET_NULL)
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