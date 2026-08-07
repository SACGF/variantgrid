"""
    Cohort-keyed annotation-version stats. The genotype-level counterpart lives
    in snpdb/models/models_cohort_stats.py to keep dependency direction
    (annotation depends on snpdb, never the other way around).
"""
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from django_extensions.db.models import TimeStampedModel

from annotation.models.models import VariantAnnotationVersion, ClinVarVersion, GeneAnnotationVersion
from snpdb.models import SampleStatsCodeVersion


class CohortGenotypeVariantAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="variant_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)

    variant_dbsnp_count = models.IntegerField(default=0)
    insertions_dbsnp_count = models.IntegerField(default=0)
    snp_dbsnp_count = models.IntegerField(default=0)
    deletions_dbsnp_count = models.IntegerField(default=0)
    ref_high_or_moderate_count = models.IntegerField(default=0)
    het_high_or_moderate_count = models.IntegerField(default=0)
    hom_high_or_moderate_count = models.IntegerField(default=0)
    unk_high_or_moderate_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "variant_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_var_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "variant_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_var_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "variant_annotation_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_var_ann_stats_aggregate_filter_uniq",
            ),
        ]

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        # Mirrors AbstractVariantAnnotationStats.count_for_zygosity — only IMPACT_HIGH_OR_MODERATE
        # is ever asked for here.
        total = 0
        if zygosity_ref:
            total += self.ref_high_or_moderate_count
        if zygosity_het:
            total += self.het_high_or_moderate_count
        if zygosity_hom:
            total += self.hom_high_or_moderate_count
        if zygosity_unk:
            total += self.unk_high_or_moderate_count
        return total


class CohortGenotypeGeneAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="gene_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    gene_annotation_version = models.ForeignKey(GeneAnnotationVersion, on_delete=CASCADE)

    gene_count = models.IntegerField(default=0)
    ref_omim_phenotype_count = models.IntegerField(default=0)
    het_omim_phenotype_count = models.IntegerField(default=0)
    hom_omim_phenotype_count = models.IntegerField(default=0)
    unk_omim_phenotype_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "gene_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_gene_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "gene_annotation_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_gene_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "gene_annotation_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_gene_ann_stats_aggregate_filter_uniq",
            ),
        ]

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        total = 0
        if zygosity_ref:
            total += self.ref_omim_phenotype_count
        if zygosity_het:
            total += self.het_omim_phenotype_count
        if zygosity_hom:
            total += self.hom_omim_phenotype_count
        if zygosity_unk:
            total += self.unk_omim_phenotype_count
        return total


class CohortGenotypeClinVarAnnotationStats(TimeStampedModel):
    cohort_genotype_collection = models.ForeignKey(
        "snpdb.CohortGenotypeCollection", on_delete=CASCADE, related_name="clinvar_annotation_stats")
    sample = models.ForeignKey("snpdb.Sample", null=True, on_delete=CASCADE)
    filter_key = models.TextField(null=True)
    passing_filter = models.BooleanField(default=False)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    clinvar_version = models.ForeignKey(ClinVarVersion, on_delete=CASCADE)

    clinvar_count = models.IntegerField(default=0)
    ref_clinvar_pathogenic_count = models.IntegerField(default=0)
    het_clinvar_pathogenic_count = models.IntegerField(default=0)
    hom_clinvar_pathogenic_count = models.IntegerField(default=0)
    unk_clinvar_pathogenic_count = models.IntegerField(default=0)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "sample", "clinvar_version", "passing_filter"],
                condition=Q(sample__isnull=False) & Q(filter_key__isnull=True),
                name="cohort_genotype_clinvar_ann_stats_per_sample_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "clinvar_version", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=True),
                name="cohort_genotype_clinvar_ann_stats_aggregate_uniq",
            ),
            models.UniqueConstraint(
                fields=["cohort_genotype_collection", "clinvar_version", "filter_key", "passing_filter"],
                condition=Q(sample__isnull=True) & Q(filter_key__isnull=False),
                name="cohort_genotype_clinvar_ann_stats_aggregate_filter_uniq",
            ),
        ]

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        total = 0
        if zygosity_ref:
            total += self.ref_clinvar_pathogenic_count
        if zygosity_het:
            total += self.het_clinvar_pathogenic_count
        if zygosity_hom:
            total += self.hom_clinvar_pathogenic_count
        if zygosity_unk:
            total += self.unk_clinvar_pathogenic_count
        return total
