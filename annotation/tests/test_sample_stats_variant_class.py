"""
Sample stats classify each call by variant class. A gene fusion's alt is symbolic but is neither an
insertion nor a deletion, and it is a class of its own rather than the SNP fallback.

A sample's totals count only its own variant calls - in a multi-sample VCF a row it is hom ref, no-call
or missing for is another sample's variant, which still shows in its per-zygosity tallies (#2087).
"""
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from annotation.tasks.calculate_sample_stats import calculate_cohort_stats
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from patients.models_enums import Zygosity
from snpdb.fake_data import create_fake_cohort
from snpdb.models import (
    VCF,
    CohortGenotype,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class SampleStatsVariantClassTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='stats_class_test')[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.sample = cls.cohort.cohortsample_set.get(sample__name="proband").sample

    def _add_genotype(self, variant, samples_zygosity: str):
        """ samples_zygosity is proband, mother, father """
        cgc = CohortGenotypeCollection.objects.get(cohort=self.cohort)
        n = len(samples_zygosity)
        CohortGenotype.objects.create(collection=cgc, variant=variant,
                                      het_count=samples_zygosity.count(Zygosity.HET),
                                      samples_zygosity=samples_zygosity,
                                      samples_allele_depth=[20] * n, samples_allele_frequency=[50] * n,
                                      samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
                                      samples_phred_likelihood=[0] * n)

    def _stats_by_sample_name(self) -> dict[str, CohortGenotypeStats]:
        calculate_cohort_stats(self.cohort, self.annotation_version)
        stats_qs = CohortGenotypeStats.objects.filter(cohort_genotype_collection__cohort=self.cohort,
                                                      sample__isnull=False, passing_filter=False)
        return {stats.sample.name: stats for stats in stats_qs.select_related("sample")}

    def test_fusion_counted_as_fusion_not_snp(self):
        proband_het = Zygosity.HET + Zygosity.MISSING + Zygosity.MISSING
        self._add_genotype(slowly_create_test_variant("1", 1000, "A", "T", self.genome_build), proband_het)
        self._add_genotype(create_gene_fusion("BCR", "ABL1").variant, proband_het)

        stats = self._stats_by_sample_name()["proband"]
        self.assertEqual(2, stats.variant_count)
        self.assertEqual(1, stats.snp_count)
        self.assertEqual(1, stats.fusions_count)

    def test_totals_only_count_the_samples_own_calls(self):
        x_snp = slowly_create_test_variant("X", 1000, "A", "T", self.genome_build)
        insertion = slowly_create_test_variant("1", 2000, "A", "AT", self.genome_build)
        self._add_genotype(x_snp, Zygosity.HET + Zygosity.HOM_REF + Zygosity.UNKNOWN_ZYGOSITY)
        self._add_genotype(insertion, Zygosity.HOM_ALT + Zygosity.HET + Zygosity.MISSING)

        stats = self._stats_by_sample_name()
        proband, mother, father = stats["proband"], stats["mother"], stats["father"]
        self.assertEqual((2, 1, 1), (proband.variant_count, proband.snp_count, proband.insertions_count))
        self.assertEqual((1, 0, 1), (mother.variant_count, mother.snp_count, mother.insertions_count))
        self.assertEqual((0, 0, 0), (father.variant_count, father.snp_count, father.insertions_count))

        self.assertEqual((1, 1), (mother.ref_count, mother.het_count))
        self.assertEqual(2, father.unk_count)
        # A hom ref call isn't chrX hom - it would raise a joint-called female's hom/het ratio
        self.assertEqual((0, 0), (mother.x_hom_count, mother.x_het_count))
        self.assertEqual(1, proband.x_het_count)

    def test_no_genotype_field_counts_unknown_zygosity_calls(self):
        VCF.objects.filter(cohort=self.cohort).update(genotype_field=None)
        no_gt_call = Zygosity.UNKNOWN_ZYGOSITY * 3
        self._add_genotype(slowly_create_test_variant("1", 3000, "C", "G", self.genome_build), no_gt_call)

        stats = self._stats_by_sample_name()["proband"]
        self.assertEqual(1, stats.variant_count)
        self.assertEqual(1, stats.snp_count)
