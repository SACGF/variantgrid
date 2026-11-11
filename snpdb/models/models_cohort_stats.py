"""
    Cohort-keyed stats models that replace the per-Sample stats family.

    Each row is keyed by CohortGenotypeCollection + nullable sample + nullable filter_key:
      sample IS NOT NULL → per-sample row (FK'd to the VCF cohort's CGC)
      sample IS NULL     → aggregate row across the cohort
      filter_key IS NULL → no extra filter beyond cohort/sample identity
      filter_key IS NOT NULL → canonical-JSON string identifying a node-side filter
                               config the calc function pre-populated (e.g. trio inheritance)
    A sub cohort has no CGC of its own: its per-sample rows are its parent's, and its aggregate rows sit
    on the parent's CGC under a filter_key carrying its samples (Cohort.get_stats_filter_key).

    See annotation/models/models_cohort_stats.py for the annotation-version-keyed
    counterparts.
"""
from typing import Optional

from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from model_utils.models import TimeStampedModel

from patients.enums import Sex
from snpdb.enums import ImportStatus, VariantContent
from snpdb.models.models_vcf import SampleStatsCodeVersion

# Targeted panels can have too few chrX calls for the hom/het ratio to mean anything
MIN_CHRX_VARIANTS_FOR_SEX_GUESS = 30
# A female's chrX has about as many hom alt calls as het, a male's hets are only pseudoautosomal or
# miscalled so hom outnumbers them many times over. The gap between is unknown: we only call when sure
CHRX_FEMALE_MAX_HOM_HET_RATIO = 1.5
CHRX_MALE_MIN_HOM_HET_RATIO = 5

DETECTED_SEX_HELP = (
    "Sex called from the sample's chrX hom/het genotype ratio - males are hemizygous so have almost "
    f"no heterozygous chrX calls. Above {CHRX_MALE_MIN_HOM_HET_RATIO} is male, below "
    f"{CHRX_FEMALE_MAX_HOM_HET_RATIO} is female, in between is unknown, and it takes at least "
    f"{MIN_CHRX_VARIANTS_FOR_SEX_GUESS} chrX calls to call it at all. Structural variants aren't counted."
)


class CohortGenotypeStats(TimeStampedModel):
    """ Genotype-level counts derived from packed CohortGenotype data —
        zygosity, variant class, x_het/x_hom. """
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
    fusions_count = models.IntegerField(default=0)  # gene-level events, @see snpdb.gene_level_variants
    # The class counts above fold these in with SNVs and indels, so they're counted again on their own
    structural_count = models.IntegerField(default=0)  # symbolic alts with coordinates (<DEL>, <DUP>, <INV>...)
    gene_copy_number_count = models.IntegerField(default=0)
    gene_splice_count = models.IntegerField(default=0)
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

    # The first SAMPLE_STATS_CODE_VERSION to count structural_count / gene_copy_number_count / gene_splice_count
    VARIANT_CONTENTS_STATS_VERSION = 8

    @property
    def total_count(self):
        return self.variant_count

    def get_variant_contents(self, has_copy_number: bool) -> dict[VariantContent, int]:
        """ Count per kind of call, leaving out the kinds with none. has_copy_number (the VCF declares a
            copy number field) is what makes its structural variants a CNV caller's rather than an SV caller's """
        gene_level_count = self.fusions_count + self.gene_copy_number_count + self.gene_splice_count
        counts = {
            VariantContent.SMALL: self.variant_count - self.structural_count - gene_level_count,
            VariantContent.CNV if has_copy_number else VariantContent.SV: self.structural_count,
            VariantContent.FUSION: self.fusions_count,
            VariantContent.GENE_COPY_NUMBER: self.gene_copy_number_count,
            VariantContent.SPLICE: self.gene_splice_count,
        }
        return {content: count for content, count in counts.items() if count}

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
    def chrx_genotyped_count(self) -> int:
        """ chrX calls the sex guess is made on - unknown-zygosity calls tell us nothing """
        return self.x_hom_count + self.x_het_count

    @property
    def chrx_hom_het_ratio(self) -> Optional[float]:
        if self.x_het_count:
            return self.x_hom_count / self.x_het_count
        return None

    @property
    def chrx_sex_guess(self) -> Sex:
        """ Genetic sex from the chrX hom/het ratio - UNKNOWN if we don't have enough calls to tell """
        if self.x_het_count and self.x_hom_count:
            if self.chrx_genotyped_count >= MIN_CHRX_VARIANTS_FOR_SEX_GUESS:
                ratio = self.chrx_hom_het_ratio
                if ratio < CHRX_FEMALE_MAX_HOM_HET_RATIO:
                    return Sex.FEMALE
                if ratio > CHRX_MALE_MIN_HOM_HET_RATIO:
                    return Sex.MALE
        return Sex.UNKNOWN

    @property
    def chrx_sex_detail(self) -> str:
        """ The chrX counts chrx_sex_guess was made on - shown alongside DETECTED_SEX_HELP """
        detail = f"This sample: {self.x_hom_count} hom / {self.x_het_count} het"
        ratio = self.chrx_hom_het_ratio
        if ratio is not None:
            detail += f", ratio {ratio:.2f}"
        return detail
