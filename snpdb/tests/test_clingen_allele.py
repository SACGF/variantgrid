from django.test import TestCase

from snpdb.clingen_allele import ClinGenAlleleServerException, get_clingen_allele, \
    ClinGenAlleleAPIException, get_clingen_allele_for_variant, variant_allele_clingen
from snpdb.models import Sequence, GenomeBuild, Variant, Locus
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI, MockServerErrorClinGenAlleleRegistryAPI
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class ClinGenAlleleTestCase(TestCase):

    def test_server_exception(self):
        """ Server responds with error code """

        with self.assertRaises(ClinGenAlleleServerException):
            get_clingen_allele("CA000072", clingen_api=MockServerErrorClinGenAlleleRegistryAPI())

    def test_api_exception(self):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        # rs6025 is reference base in 38, so gives - NoConsistentAlignment
        variant_rs6025 = slowly_create_test_variant("1", 169519049, "T", "C", grch37)

        with self.assertRaises(ClinGenAlleleAPIException):
            get_clingen_allele_for_variant(grch37, variant_rs6025, clingen_api=MockClinGenAlleleRegistryAPI())

    def _get_variant_allele_failed_clingen(self):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        variant = slowly_create_test_variant("3", 128198980, 'A', 'T', grch37)

        clingen_api_fail = MockServerErrorClinGenAlleleRegistryAPI()
        return variant_allele_clingen(grch37, variant, clingen_api=clingen_api_fail)

    def test_fail_then_retry_success(self):
        variant_allele = self._get_variant_allele_failed_clingen()
        self.assertIsNone(variant_allele.allele.clingen_allele, "Allele.clingen_allele not set after API failure")
        # Attempt to retry this time with success
        clingen_api_success = MockClinGenAlleleRegistryAPI()
        variant_allele = variant_allele_clingen(variant_allele.genome_build, variant_allele.variant,
                                                existing_variant_allele=variant_allele,
                                                clingen_api=clingen_api_success)
        self.assertIsNotNone(variant_allele.allele.clingen_allele, "Allele.clingen_allele set after API success")
        self.assertIsNone(variant_allele.error, "VariantAllele.error cleared after API success")

    def test_fail_then_retry_success_existing_allele(self):
        # Used to have a "ValueError: Cannot assign" due to assigning an Allele field a VariantAllele
        clingen_api_success = MockClinGenAlleleRegistryAPI()
        # Existing Allele w/clingen_allele - should be merged
        clingen_allele = get_clingen_allele("CA10617208", clingen_api=clingen_api_success)
        variant_allele = self._get_variant_allele_failed_clingen()

        variant_allele = variant_allele_clingen(variant_allele.genome_build, variant_allele.variant,
                                                existing_variant_allele=variant_allele,
                                                clingen_api=clingen_api_success)
        self.assertEqual(clingen_allele.allele, variant_allele.allele, "Alleles merged")
