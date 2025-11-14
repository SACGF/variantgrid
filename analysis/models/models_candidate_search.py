"""

Find/surface potential variants / classifications for human review

Note on naming:
* Review, Flag are already used in VG code base
* "Finding" in clinical contexts means real observed/reported result

"""
from django.contrib.auth.models import User
from django.db import models
from django.db.models import CASCADE
from django.db.models.deletion import SET_NULL
from model_utils.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from analysis.models import Analysis
from annotation.models import ClinVar, AnnotationVersion
from classification.models import Classification
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from patients.models_enums import Zygosity
from snpdb.models import Variant, Sample


class CandidateSearchType(models.TextChoices):
    REANALYSIS_NEW_ANNOTATION = 'R', 'Reanalysis - new annotation'
    CROSS_SAMPLE_CLASSIFICATION = 'S', 'Cross Sample Classification'
    CLASSIFICATION_EVIDENCE_UPDATE = 'E', 'Classification Evidence Update'


class CandidateStatus(models.TextChoices):
    OPEN = 'O', 'OPEN'
    RESOLVED = 'R', 'RESOLVED'
    HIDDEN = 'D', 'HIDDEN'
    HIGHLIGHTED = 'L', 'HIGHLIGHTED'


class CandidateSearchVersion(TimeStampedModel):
    search_type = models.CharField(choices=CandidateSearchType.choices, max_length=1)
    code_version = models.IntegerField()
    class Meta:
        unique_together = ('search_type', 'code_version')


class CandidateSearchRun(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    search_version = models.ForeignKey(CandidateSearchVersion, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    config_snapshot = models.JSONField(default=dict)
    git_hash = models.TextField()


class AbstractCandidate(TimeStampedModel):
    objects = InheritanceManager()
    search_run = models.ForeignKey(CandidateSearchRun, on_delete=CASCADE)
    status = models.CharField(choices=CandidateStatus.choices, max_length=1, default=CandidateStatus.OPEN)
    notes = models.TextField(null=True, blank=True)
    evidence = models.JSONField(default=dict)
    reviewer = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    reviewer_comment = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True
        ordering = ('-created',)


class ReanalysisCandidate(AbstractCandidate):
    """ Surface a variant in an analysis, using newer annotation """
    # Analysis. Annotation used is analysis.annotation_version
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    # Annotation used to find updates
    annotation_version = models.ForeignKey(AnnotationVersion, on_delete=CASCADE)
    clinvar = models.ForeignKey(ClinVar, null=True, blank=True, on_delete=CASCADE)
    # Not currently used
    variant = models.ForeignKey(Variant, null=True, blank=True, on_delete=CASCADE)


class AbstractClassificationCandidate(AbstractCandidate):
    classification = models.ForeignKey(Classification, on_delete=CASCADE)

    class Meta:
        abstract = True
        ordering = ('-created',)

class CrossSampleClassificationCandidate(AbstractClassificationCandidate):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    zygosity = models.CharField(choices=Zygosity.CHOICES, max_length=1)


class ClassificationEvidenceUpdateCandidate(AbstractClassificationCandidate):
    """ Examples:
            * ClinVar appears for classification
            * Splice AI calculated for a VUS
    """
    annotation_version = models.ForeignKey(AnnotationVersion, null=True, on_delete=CASCADE)
