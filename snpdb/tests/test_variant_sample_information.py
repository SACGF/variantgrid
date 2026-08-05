from django.contrib.auth.models import User
from django.test import TestCase

from patients.models_enums import Zygosity
from snpdb.models import CohortGenotype, CohortGenotypeCollection, GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variant_sample_information import VariantZygosityCounts


class VariantZygosityCountsTest(TestCase):
    """ Sample scoping comes from the variant's contig, which can belong to more than one build """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='vsi_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

        # A cohort/VCF per build - each has a proband the user can read
        cls.cohort_37 = create_fake_cohort(cls.user, cls.grch37)
        cls.cohort_38 = create_fake_cohort(cls.user, cls.grch38)

    def _het_proband(self, cohort, variant):
        """ samples_zygosity is indexed [proband, mother, father] - only proband is called """
        cgc = CohortGenotypeCollection.objects.get(cohort=cohort)
        samples_zygosity = Zygosity.HET + Zygosity.MISSING + Zygosity.MISSING
        n = len(samples_zygosity)
        CohortGenotype.objects.create(collection=cgc, variant=variant, het_count=1,
                                      samples_zygosity=samples_zygosity,
                                      samples_allele_depth=[20] * n, samples_allele_frequency=[100] * n,
                                      samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
                                      samples_phred_likelihood=[0] * n)

    def test_shared_contig_counts_both_builds(self):
        """ GRCh37/38 share the MT contig, so both builds' VCFs produce the same variant. Samples from
            both must be visible rows, not counted as observations the user isn't allowed to see """
        variant = slowly_create_test_variant("MT", 1000, "A", "T", self.grch37)
        self.assertEqual({self.grch37, self.grch38}, variant.genome_builds,
                         "MT contig is shared between GRCh37 and GRCh38")

        self._het_proband(self.cohort_37, variant)
        self._het_proband(self.cohort_38, variant)

        counts = VariantZygosityCounts(self.user, variant)
        self.assertEqual(2, counts.num_observations)
        self.assertEqual(2, counts.num_visible_observations)
        self.assertEqual(0, counts.num_invisible_observations)
        self.assertEqual(2, counts.num_visible_het)

    def test_build_specific_contig_only_counts_own_build(self):
        """ Normal contigs aren't shared - a GRCh38 sample is a different variant, so must not be counted """
        variant_37 = slowly_create_test_variant("3", 1000, "A", "T", self.grch37)
        variant_38 = slowly_create_test_variant("3", 1000, "A", "T", self.grch38)
        self.assertNotEqual(variant_37.pk, variant_38.pk)
        self.assertEqual({self.grch37}, variant_37.genome_builds)

        self._het_proband(self.cohort_37, variant_37)
        self._het_proband(self.cohort_38, variant_38)

        counts = VariantZygosityCounts(self.user, variant_37)
        self.assertEqual(1, counts.num_observations)
        self.assertEqual(1, counts.num_visible_observations)
        self.assertEqual(0, counts.num_invisible_observations)

    def test_hidden_samples(self):
        """ has_hidden_samples is what the template branches on to explain the visible/total split """
        variant = slowly_create_test_variant("3", 2000, "A", "T", self.grch37)

        # Owner can read the whole VCF, so nothing is hidden from them
        counts = VariantZygosityCounts(self.user, variant)
        self.assertEqual(counts.num_samples, counts.num_user_samples)
        self.assertFalse(counts.has_hidden_samples)

        other_user = User.objects.create(username='vsi_test_no_perms')
        no_perms = VariantZygosityCounts(other_user, variant)
        self.assertEqual(0, no_perms.num_user_samples)
        self.assertTrue(no_perms.has_hidden_samples)
