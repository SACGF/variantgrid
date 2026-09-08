"""
Sample stats classify each call by variant class. A gene fusion's alt is symbolic but is neither an
insertion nor a deletion, and it is a class of its own rather than the SNP fallback.
"""
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tasks.calculate_sample_stats import calculate_cohort_stats
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from patients.models_enums import Zygosity
from snpdb.models import CohortGenotype, CohortGenotypeCollection, CohortGenotypeStats, GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
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

    def _proband_het(self, variant):
        cgc = CohortGenotypeCollection.objects.get(cohort=self.cohort)
        samples_zygosity = Zygosity.HET + Zygosity.MISSING + Zygosity.MISSING
        n = len(samples_zygosity)
        CohortGenotype.objects.create(collection=cgc, variant=variant, het_count=1,
                                      samples_zygosity=samples_zygosity,
                                      samples_allele_depth=[20] * n, samples_allele_frequency=[50] * n,
                                      samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
                                      samples_phred_likelihood=[0] * n)

    def test_fusion_counted_as_fusion_not_snp(self):
        self._proband_het(slowly_create_test_variant("1", 1000, "A", "T", self.genome_build))
        self._proband_het(create_gene_fusion("BCR", "ABL1").variant)

        calculate_cohort_stats(self.cohort, self.annotation_version)

        stats = CohortGenotypeStats.objects.get(sample=self.sample, passing_filter=False)
        self.assertEqual(2, stats.variant_count)
        self.assertEqual(1, stats.snp_count)
        self.assertEqual(1, stats.fusions_count)
