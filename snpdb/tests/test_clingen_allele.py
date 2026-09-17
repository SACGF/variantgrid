from unittest.mock import patch

from django.db import IntegrityError
from django.test import TestCase, override_settings

from annotation.fake_annotation import get_fake_annotation_version
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.clingen_allele import (
    ClinGenAlleleAPIException,
    ClinGenAlleleServerException,
    _create_variant_allele_with_new_allele,
    get_clingen_allele,
    get_clingen_allele_for_variant,
    get_variant_allele_for_variant,
    populate_clingen_alleles_for_variants,
    variant_allele_clingen,
)
from snpdb.models import Allele, ClinGenAllele, GenomeBuild, VariantAllele, VariantCoordinate
from snpdb.tests.utils.mock_clingen_api import (
    MockClinGenAlleleRegistryAPI,
    MockServerErrorClinGenAlleleRegistryAPI,
)
from snpdb.tests.utils.vcf_testing_utils import (
    create_mock_allele,
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)


class ClinGenAlleleTestCase(TestCase):

    @classmethod
    def setUpTestData(cls):
        # We need this as HGVSMatcher (biocommons) needs annotation version / genomes
        super().setUpTestData()
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)

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
        self.assertIsNone(variant_allele.clingen_error, "VariantAllele.clingen_error cleared after API success")

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


class ClinGenAlleleNeverRegisteredTestCase(TestCase):
    """ Variants that can never get a ClinGenAllele used to collect a new Allele per call - #1361 / #1844 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        # 1bp over the ClinGen limit, so clingen_allele_skip_reason() will always refuse it
        variant_coordinate = VariantCoordinate(chrom="3", position=128198980, ref="A",
                                               alt=VCFSymbolicAllele.DEL,
                                               svlen=-(ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE + 1))
        cls.variant = slowly_create_test_variant_from_coordinate(variant_coordinate, cls.grch37)

    def _populate(self):
        populate_clingen_alleles_for_variants(self.grch37, [self.variant],
                                              clingen_api=MockClinGenAlleleRegistryAPI())

    def test_populate_creates_one_allele_for_new_never_clingen_variant(self):
        self.assertIsNotNone(self.variant.clingen_allele_skip_reason())
        self._populate()
        variant_allele = VariantAllele.objects.get(variant=self.variant, genome_build=self.grch37)
        self.assertIsNone(variant_allele.clingen_error)
        self.assertIsNone(variant_allele.allele.clingen_allele)

    def test_populate_never_clingen_variant_is_idempotent(self):
        for _ in range(3):
            self._populate()

        va_qs = VariantAllele.objects.filter(variant=self.variant, genome_build=self.grch37)
        self.assertEqual(va_qs.count(), 1)
        self.assertEqual(Allele.objects.filter(variantallele__variant=self.variant).distinct().count(), 1)

    def test_populate_repeated_variant_in_one_call(self):
        populate_clingen_alleles_for_variants(self.grch37, [self.variant, self.variant],
                                              clingen_api=MockClinGenAlleleRegistryAPI())
        self.assertEqual(VariantAllele.objects.filter(variant=self.variant).count(), 1)
        self.assertEqual(Allele.objects.count(), 1)

    @override_settings(CLINGEN_ALLELE_REGISTRY_LOGIN=None)
    def test_get_variant_allele_for_variant_creates_then_reuses(self):
        variant_allele = get_variant_allele_for_variant(self.grch37, self.variant)
        self.assertEqual(Allele.objects.count(), 1)

        again = get_variant_allele_for_variant(self.grch37, self.variant)
        self.assertEqual(again.pk, variant_allele.pk)
        self.assertEqual(Allele.objects.count(), 1)

    @override_settings(CLINGEN_ALLELE_REGISTRY_LOGIN=None)
    def test_create_variant_allele_uses_the_race_winners_link(self):
        winning_allele = create_mock_allele(self.variant, self.grch37)
        with patch.object(VariantAllele.objects, "create", side_effect=IntegrityError):
            variant_allele = _create_variant_allele_with_new_allele(self.variant, self.grch37)

        self.assertEqual(variant_allele.allele, winning_allele)
        # The Allele we made before losing the race is gone, rather than left with nothing pointing at it
        self.assertEqual(Allele.objects.count(), 1)
