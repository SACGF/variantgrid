from django.db.models import IntegerChoices, TextChoices
from django.db import models
from django_extensions.db.models import TimeStampedModel

from classification.enums import OverlapStatus, TestingContextBucket
from classification.models import ClassificationGrouping
from snpdb.models import Allele, Lab


class OverlapType(TextChoices):
    SINGLE_CONTEXT = "context", "Single Context"
    CROSS_CONTEXT = "cross", "Cross Context"
    CLINVAR_EXPERT_PANEL = "clinvar", "ClinVar Expert Panel"


class ClassificationResultValue(TextChoices):
    ONC_PATH = "O", "Onco-Path"
    CLINICAL_SIGNIFICANCE = "S", "Clinical significance"


class OverlapContributionStatus(TextChoices):
    PENDING_CALCULATION = "P", "Pending Calculation"
    CONTRIBUTING = "C", "Contributing"
    NOT_SHARED = "N", "Not-shared"
    NO_VALUE = "X", "No value"
    NON_COMPARABLE_VALUE = "Z", "Non-comparable value"  # e.g. Risk Factor


class TriageStatus(TextChoices):
    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Classification"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"


class Overlap(TimeStampedModel):
    """
    Overlap is made by composition as making a separate model for each overlap type added a lot of overhead for just issolating a few fields
    """
    overlap_type = models.TextField(choices=OverlapType.choices)
    value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
    testing_context = models.TextField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
    tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
    lab = models.ForeignKey(Lab, on_delete=models.CASCADE, null=True, blank=True)  # only use for

    overlap_status = models.TextField(max_length=1, choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)

    class Meta:
        indexes = [models.Index(fields=['overlap_type']), models.Index(fields=['allele'])]
        unique_together = ('overlap_type', 'allele', 'value_type', 'testing_context', 'tumor_type_category', 'lab')


class ClassificationGroupingOverlapContribution(TimeStampedModel):
    classification_grouping = models.ForeignKey(ClassificationGrouping, on_delete=models.CASCADE)
    overlap = models.ForeignKey(Overlap, on_delete=models.CASCADE)
    contribution_status = models.TextField(max_length=1, choices=OverlapContributionStatus.choices, default=OverlapContributionStatus.PENDING_CALCULATION)

    class Meta:
        unique_together = ('classification_grouping', 'overlap')


class ClassificationGroupingValueTriage(TimeStampedModel):
    classification_grouping = models.ForeignKey(ClassificationGrouping, on_delete=models.CASCADE)
    result_value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    new_value = models.TextField(null=True, blank=True)
    triage_status = models.TextField(max_length=1, choices=TriageStatus.choices, default=TriageStatus.PENDING)

    class Meta:
        unique_together = ('classification_grouping', 'result_value_type')