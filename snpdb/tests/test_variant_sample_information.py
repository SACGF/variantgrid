from django.contrib.auth.models import User
from django.test import TestCase

from patients.models_enums import Zygosity
from snpdb.models import CohortGenotype, CohortGenotypeCollection, GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variant_sample_information import VariantSampleGenotypes, VariantZygosityCounts


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
        self._proband_zygosity(cohort, variant, Zygosity.HET)

    def _proband_zygosity(self, cohort, variant, zygosity):
        cgc = CohortGenotypeCollection.objects.get(cohort=cohort)
        samples_zygosity = zygosity + Zygosity.MISSING + Zygosity.MISSING
        n = len(samples_zygosity)
        CohortGenotype.objects.create(collection=cgc, variant=variant,
                                      het_count=int(zygosity == Zygosity.HET),
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


class LocusCountsTest(TestCase):
    """ The locus counts table is per-allele zygosity counts for every variant at the locus """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='locus_counts_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)

    def _proband_zygosity(self, variant, zygosity):
        cgc = CohortGenotypeCollection.objects.get(cohort=self.cohort)
        samples_zygosity = zygosity + Zygosity.MISSING + Zygosity.MISSING
        n = len(samples_zygosity)
        CohortGenotype.objects.create(collection=cgc, variant=variant,
                                      het_count=int(zygosity == Zygosity.HET),
                                      samples_zygosity=samples_zygosity,
                                      samples_allele_depth=[20] * n, samples_allele_frequency=[100] * n,
                                      samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
                                      samples_phred_likelihood=[0] * n)

    def _locus_counts_by_variant(self, variant) -> dict:
        locus_counts = VariantSampleGenotypes(self.user, variant).to_json()["locus_counts"]
        return {row["variant"]: row for row in locus_counts}

    def test_allele_without_samples_counts_zero(self):
        """ An allele at the locus with no CohortGenotype still gets a row, but it has no observations """
        variant = slowly_create_test_variant("3", 3000, "A", "T", self.grch37)
        no_samples = slowly_create_test_variant("3", 3000, "A", "G", self.grch37)
        self._proband_zygosity(variant, Zygosity.HET)

        by_variant = self._locus_counts_by_variant(variant)
        self.assertEqual({str(variant), str(no_samples)}, set(by_variant))

        row = by_variant[str(no_samples)]
        self.assertEqual(0, row["total"])
        self.assertEqual(0, row["Unknown"], "No CohortGenotype isn't an unknown zygosity call")

    def test_unknown_zygosity_is_counted(self):
        """ './.' calls are real observations, and the grid has a zygosity filter for them """
        variant = slowly_create_test_variant("3", 4000, "A", "T", self.grch37)
        self._proband_zygosity(variant, Zygosity.UNKNOWN_ZYGOSITY)

        data = VariantSampleGenotypes(self.user, variant).to_json()
        row = {r["variant"]: r for r in data["locus_counts"]}[str(variant)]
        self.assertEqual(1, row["total"])
        self.assertEqual(1, row["Unknown"])
        self.assertEqual({Zygosity.UNKNOWN_ZYGOSITY: 1}, data["zygosity_counts"])
