from django.db import models
from django.db.models.deletion import CASCADE
from model_utils.managers import InheritanceManager

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models.models import VariantAnnotationVersion, VariantAnnotation
from genes.models import Gene
from library.utils import rgb_invert
from snpdb.models import Sample, Cohort, ShareLevel, GenomeBuild, Variant
from snpdb.models.models_enums import ProcessingStatus
from classification.models import Classification, ClassificationModification
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

    @staticmethod
    def handle_classification_change(classification):
        sample = classification.sample
        if sample:
            gene_by_vav_id = GeneCountType.get_gene_by_variant_annotation_version(classification)
            for gene_count_type in GeneCountType.objects.filter(enabled=True, uses_classifications=True):
                # Search through all modifications (not just last published) - so you recount if ever counted
                qs = gene_count_type.get_classification_qs(last_published=False)
                if qs.filter(pk=classification.pk).exists():
                    print(f"Need to recalculate {gene_count_type} - for {classification}")
                    for cgc in gene_count_type.cohortgenecounts_set.filter(cohort__cohortsample__sample=sample):
                        gene = gene_by_vav_id.get(cgc.variant_annotation_version.pk)
                        if gene:
                            cgc.launch_task(gene.pk)

    @staticmethod
    def get_gene_by_variant_annotation_version(classification):
        gene_by_vav_id = {}
        for vav in VariantAnnotationVersion.objects.all():
            variant_annotation = classification.get_variant_annotation(vav)
            if variant_annotation:
                gene_by_vav_id[vav.pk] = variant_annotation.gene
        return gene_by_vav_id

    def _get_variant_q(self, genome_build: GenomeBuild):
        q_variant = None
        if self.uses_classifications:
            vc_qs = self.get_classification_qs()
            q_variant = Classification.get_variant_q_from_classification_qs(vc_qs, genome_build)
        return q_variant

    def get_variant_queryset(self, variant_annotation_version: VariantAnnotationVersion):
        annotation_version = variant_annotation_version.get_any_annotation_version()

        qs = get_variant_queryset_for_annotation_version(annotation_version=annotation_version)
        kwargs = {VariantAnnotation.GENE_COLUMN + "__isnull": False}
        qs = qs.filter(Variant.get_no_reference_q(), **kwargs)
        q_variant = self._get_variant_q(annotation_version.genome_build)
        if q_variant:
            qs = qs.filter(q_variant)
        return qs

    def get_classification_qs(self, last_published: bool = True):
        """ @param last_published - set to False to look through historical modifications """

        if not self.uses_classifications:
            raise ValueError(f"{self} does has uses_classifications=False")

        if self.name == "RUNX1_classified_damage":
            kwargs = {"share_level__in": ShareLevel.same_and_higher(ShareLevel.ALL_USERS),
                      "published_evidence__clinical_significance__note__icontains": 'Clinically Relevant',
                      "published_evidence__allele_origin__value__in": ['somatic', 'likely_somatic'],
                      "published_evidence__specimen_id__value__icontains": 'ML'}
            if last_published:
                kwargs["is_last_published"] = True
            vcm_qs = ClassificationModification.objects.filter(**kwargs)
            qs = Classification.objects.filter(pk__in=vcm_qs.values('classification'))
        else:
            qs = Classification.objects.none()
        return qs

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
    processing_status = models.CharField(max_length=1, choices=ProcessingStatus.CHOICES,
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
