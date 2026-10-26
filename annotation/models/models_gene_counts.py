from django.db import models
from django.db.models.deletion import CASCADE
from model_utils.managers import InheritanceManager

from annotation.models.models_version import VariantAnnotationVersion
from genes.models import Gene
from library.utils import rgb_invert
from snpdb.enums import ProcessingStatus
from snpdb.models import Cohort, Sample
from variantgrid.celery import app


class VariantSource(models.Model):
    objects = InheritanceManager()

    def get_queryset(self):
        return None


class SampleAnnotationVersionVariantSource(VariantSource):
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ("sample", "variant_annotation_version")

    def __str__(self):
        return f"VariantSource from sample_id: {self.sample_id}, vav: {self.variant_annotation_version}"


class GeneCountType(models.Model):
    name = models.TextField(primary_key=True)
    celery_task_name = models.TextField(null=True)
    enabled = models.BooleanField(default=False)
    uses_classifications = models.BooleanField(default=False)

    def __str__(self):
        return self.name


class GeneValue(models.Model):
    gene_count_type = models.ForeignKey(GeneCountType, on_delete=CASCADE)
    label = models.TextField()
    rgb = models.CharField(max_length=7)  # '#rrggbb'
    show_counts = models.BooleanField(default=True)
    use_as_empty_value = models.BooleanField(default=False)  # Only 1 per gene_count_type

    class Meta:
        unique_together = ('gene_count_type', 'label')

    @property
    def inverted_rgb(self):
        return rgb_invert(self.rgb)

    @property
    def style(self):
        return f"color: {self.inverted_rgb}; background-color: {self.rgb};"

    def __str__(self):
        return f"{self.gene_count_type}: {self.label}"


class GeneValueCountCollection(models.Model):
    source = models.ForeignKey(VariantSource, on_delete=CASCADE)
    gene_count_type = models.ForeignKey(GeneCountType, on_delete=CASCADE)

    def __str__(self):
        return f"{self.source}: {self.gene_count_type}"


class GeneValueCount(models.Model):
    """ A way to count hits for a gene, so we can draw a sample_gene_matrix """
    collection = models.ForeignKey(GeneValueCountCollection, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    value = models.ForeignKey(GeneValue, on_delete=CASCADE)
    count = models.IntegerField()

    class Meta:
        unique_together = ("collection", "gene", "value")

    def __str__(self):
        return f"{self.collection}/{self.gene}/{self.value}: {self.count}"


class CohortGeneCounts(models.Model):
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    gene_count_type = models.ForeignKey(GeneCountType, on_delete=CASCADE)
    processing_status = models.CharField(max_length=1, choices=ProcessingStatus.choices,
                                         default=ProcessingStatus.CREATED)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()
    celery_task = models.CharField(max_length=36, null=True)

    def launch_task(self, *extra_args):
        args = (self.pk,) + extra_args
        result = app.send_task(self.gene_count_type.celery_task_name, args=args)
        self.celery_task = result.id
        self.save()
        return self.celery_task

    def __str__(self):
        processing_status = self.get_processing_status_display()
        return f"{self.cohort} - {self.gene_count_type} ({processing_status})"
