from typing import Optional
from django.db import models
from model_utils.models import TimeStampedModel
from classification.json_utils import ValidatedJson
from classification.models import ClassificationModification, ConditionResolved
from classification.models.clinvar_export_convertor import ClinVarExportConverter
from snpdb.models import ClinVarKey, Allele


class ClinVarAllele(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar allele"

    allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=models.CASCADE)

    def __str__(self):
        return f"{self.allele} {self.clinvar_key}"


class ClinVarExportRecord(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export record"

    clinvar_allele = models.ForeignKey(ClinVarAllele, null=True, blank=True, on_delete=models.CASCADE)
    condition = models.JSONField()
    classification_based_on = models.ForeignKey(ClassificationModification, null=True, blank=True, on_delete=models.CASCADE)
    snv = models.TextField(null=True, blank=True)  # if not set yet

    content_last_submitted = models.JSONField(null=True, blank=True)

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
        self.classification_based_on = new_classification_based_on
        # FIXME, check to see if we changed since last submission
        self.save()

    @staticmethod
    def new_condition(clinvar_allele: ClinVarAllele, condition: ConditionResolved, candidate: Optional[ClassificationModification]):
        cc = ClinVarExportRecord(clinvar_allele=clinvar_allele, condition=condition.as_json_minimal())
        cc.update_classification(candidate)

    @property
    def content_current(self) -> Optional[ValidatedJson]:
        if classification := self.classification_based_on:
            converter = ClinVarExportConverter(classification)
            return converter.as_json()
        return None


class ClinVarStatus(models.TextChoices):
    """
    For when multiple terms exist, are we uncertain which one it is, or are we curating for multiple.
    NOT_DECIDED means a choice is still required
    """
    PENDING = "P"
    SUBMITTED = "S"
    ACKNOWLEDGED = "A"
    IN_ERROR = "E"


class ClinVarExportRecordSubmission(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar export record submission"

    clinvar_candidate = models.ForeignKey(ClinVarExportRecord, on_delete=models.PROTECT)  # if there's been an actual submission, don't allow deletes
    status = models.CharField(max_length=1, choices=ClinVarStatus.choices)
    submitted_json = models.JSONField()


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