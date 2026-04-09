from django.test import TestCase

from snpdb.models import Allele, AlleleConversionTool, AlleleOrigin, GenomeBuild, VariantAllele
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


class AlleleTestCase(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_allele_merge(self):
        allele_1 = Allele.objects.create()
        allele_2 = Allele.objects.create()

        allele_1.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele_2)


class AlleleMergeEdgeCasesTestCase(TestCase):
    """
    Existing test_allele_merge only checks the happy path doesn't crash.
    These tests cover the ClinGen rejection branch and verify data integrity.
    """

    _clingen_counter = 0

    def _make_allele_with_clingen(self):
        AlleleMergeEdgeCasesTestCase._clingen_counter += 1
        n = AlleleMergeEdgeCasesTestCase._clingen_counter
        return Allele.objects.create(
            clingen_allele=ClinGenAllele.objects.create(id=n, api_response={"id": f"CA{n}"})
        )

    def test_merge_returns_false_when_both_have_clingen(self):
        allele_a = self._make_allele_with_clingen()
        allele_b = self._make_allele_with_clingen()
        self.assertFalse(allele_a.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele_b))

    def test_failed_merge_leaves_variant_alleles_on_original(self):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        allele_a = self._make_allele_with_clingen()
        allele_b = self._make_allele_with_clingen()
        v = slowly_create_test_variant("1", 5_000_000, "A", "T", grch37)
        VariantAllele.objects.create(
            variant=v, genome_build=grch37, allele=allele_b,
            origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            allele_linking_tool=AlleleConversionTool.DBSNP,
        )

        allele_a.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele_b)

        self.assertEqual(v.variantallele_set.get().allele, allele_b)

    def test_merge_self_raises(self):
        allele = Allele.objects.create()
        with self.assertRaises(ValueError):
            allele.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele)

    def test_successful_merge_moves_variant_alleles(self):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        allele_a = Allele.objects.create()
        allele_b = Allele.objects.create()
        v = slowly_create_test_variant("1", 6_000_000, "C", "G", grch37)
        VariantAllele.objects.create(
            variant=v, genome_build=grch37, allele=allele_b,
            origin=AlleleOrigin.IMPORTED_TO_DATABASE,
            allele_linking_tool=AlleleConversionTool.DBSNP,
        )

        allele_a.merge(AlleleConversionTool.BCFTOOLS_LIFTOVER, allele_b)

        self.assertEqual(v.variantallele_set.get().allele, allele_a)


class AlleleVariantForBuildTestCase(TestCase):
    """
    variant_for_build and variant_for_any_build(best_attempt=False) must raise
    ValueError when the requested build is absent. Callers depend on this contract.
    """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        v37 = slowly_create_test_variant("1", 3_000_000, "C", "G", cls.grch37)
        cls.allele = create_mock_allele(v37, cls.grch37)

    def test_variant_for_build_missing_raises(self):
        with self.assertRaises(ValueError):
            self.allele.variant_for_build(self.grch38)

    def test_variant_for_any_build_missing_best_attempt_false_raises(self):
        with self.assertRaises(ValueError):
            self.allele.variant_for_any_build(self.grch38, best_attempt=False)
