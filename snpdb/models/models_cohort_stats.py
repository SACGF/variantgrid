"""
    Cohort-keyed stats models that replace the per-Sample stats family.

    Each row is keyed by CohortGenotypeCollection + nullable sample + nullable filter_key:
      sample IS NOT NULL → per-sample row (FK'd to the VCF cohort's CGC)
      sample IS NULL     → aggregate row across the cohort
      filter_key IS NULL → no extra filter beyond cohort/sample identity
      filter_key IS NOT NULL → canonical-JSON string identifying a node-side filter
                               config the calc function pre-populated (e.g. trio inheritance)

    See annotation/models/models_cohort_stats.py for the annotation-version-keyed
    counterparts.
"""
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from django_extensions.db.models import TimeStampedModel

from patients.models_enums import Sex
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_vcf import SampleStatsCodeVersion


class CohortGenotypeStats(TimeStampedModel):
    """ Replaces SampleStats(+PassingFilter). Counts derived from packed
        CohortGenotype data — zygosity, variant class, x_het/x_hom. """
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="genotype_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices,
                                     default=ImportStatus.CREATED)

    variant_count = models.IntegerField(default=0)
    snp_count = models.IntegerField(default=0)
    insertions_count = models.IntegerField(default=0)
    deletions_count = models.IntegerField(default=0)
    ref_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)
    unk_count = models.IntegerField(default=0)
    x_hom_count = models.IntegerField(default=0)
    x_het_count = models.IntegerField(default=0)
    x_unk_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_stats_aggregate_filter_uniq",
            ),
        ]
        indexes = [models.Index(fields=["sample", "passing_filter"])]

    @property
    def total_count(self):
        return self.variant_count

    @staticmethod
    def percent(a, b):
        percent = float('NaN')
        if b:
            percent = 100.0 * a / b
        return percent

    @property
    def variant_percent(self):
        return self.percent(self.variant_count, self.total_count)

    @property
    def snp_percent(self):
        return self.percent(self.snp_count, self.total_count)

    @property
    def insertions_percent(self):
        return self.percent(self.insertions_count, self.total_count)

    @property
    def deletions_percent(self):
        return self.percent(self.deletions_count, self.total_count)

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        count = 0
        if zygosity_ref:
            count += self.ref_count
        if zygosity_het:
            count += self.het_count
        if zygosity_hom:
            count += self.hom_count
        if zygosity_unk:
            count += self.unk_count
        return count

    @property
    def chrx_sex_guess(self):
        sex = Sex.UNKNOWN
        if self.x_het_count and self.x_hom_count:
            ratio = self.x_hom_count / self.x_het_count
            if ratio < 0.2:
                sex = Sex.FEMALE
            elif ratio > 0.8:
                sex = Sex.MALE
        return Sex(sex).label
