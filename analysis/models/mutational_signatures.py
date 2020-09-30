from django.contrib.postgres.fields.array import ArrayField
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.fields import FloatField

from analysis.models.enums import MinimisationResultType, MinimisationStrategy
from snpdb.models import Sample, SoftwareVersion
from snpdb.models.models_enums import ImportStatus


class MutationalSignatureCalculator(SoftwareVersion):
    num_iterations = models.IntegerField()
    sampling_fraction = models.FloatField()
    signature_data_filename = models.TextField()
    minimisation_strategy = models.CharField(max_length=2, choices=MinimisationStrategy.CHOICES)

    class Meta:
        unique_together = ("name", "version", "num_iterations", "sampling_fraction", "signature_data_filename", "minimisation_strategy")


class MutationalSignature(models.Model):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    calculator = models.ForeignKey(MutationalSignatureCalculator, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.CHOICES, default=ImportStatus.CREATED)
    summary = models.TextField()
    mean = ArrayField(FloatField(), null=True)
    num_snps = models.IntegerField(null=True)

    class Meta:
        unique_together = ("sample", "calculator")

    def __str__(self):
        if self.import_status == ImportStatus.SUCCESS:
            description = self.summary
        else:
            import_status = self.get_import_status_display()
            description = f"MutationalSignature {self.pk}, import_status={import_status}"
        return description


class MutationalSignatureMinimisationResult(models.Model):
    mutational_signature = models.ForeignKey(MutationalSignature, on_delete=CASCADE)
    result_type = models.CharField(max_length=1, choices=MinimisationResultType.CHOICES)
    iteration = models.IntegerField()
    solution_array = ArrayField(FloatField())
    fit_data = ArrayField(FloatField())
    diff_data = ArrayField(FloatField())
    ls_sum_diff = models.FloatField()
    la_sum_diff = models.FloatField()

    class Meta:
        unique_together = ("mutational_signature", "result_type", "iteration")


class MutationalSignatureMutationCount(models.Model):
    mutational_signature = models.ForeignKey(MutationalSignature, on_delete=CASCADE)
    mutation_type = models.TextField()
    count = models.IntegerField()
