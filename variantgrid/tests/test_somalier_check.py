"""
The deployment check that keeps settings.SOMALIER["compensate_allele_order"] honest against the
somalier that is actually installed (@see snpdb/variants_to_vcf.py:somalier_alleles_flipped).
"""
from django.test import SimpleTestCase

from variantgrid.deployment_validation.somalier_check import (
    ALLELE_ORDER_SITES_PER_ORDER,
    allele_order_result,
    allele_order_unchecked,
)

NUM_SITES = 2 * ALLELE_ORDER_SITES_PER_ORDER
HALF_INVERTED = (ALLELE_ORDER_SITES_PER_ORDER, ALLELE_ORDER_SITES_PER_ORDER)


class SomalierAlleleOrderResultTest(SimpleTestCase):
    """ Every call the check writes is hom-ref, so the counts coming back say who is right """

    def test_all_hom_ref_is_valid(self):
        for compensating in (True, False):
            with self.subTest(compensating=compensating):
                self.assertTrue(allele_order_result("GRCh38", NUM_SITES, NUM_SITES, 0, compensating)["valid"])

    def test_half_inverted_says_which_way_to_set_it(self):
        """ Only the half whose alleles disagree flipping is the setting being wrong, nothing else """
        for compensating in (True, False):
            with self.subTest(compensating=compensating):
                result = allele_order_result("GRCh38", NUM_SITES, *HALF_INVERTED, compensating)
                self.assertFalse(result["valid"])
                self.assertIn(f'compensate_allele_order"] = {not compensating}', result["fix"])
                self.assertIn("somalier_existing_vcfs --clear", result["fix"])

    def test_anything_else_does_not_claim_to_know_why(self):
        result = allele_order_result("GRCh38", NUM_SITES, 0, NUM_SITES, True)
        self.assertFalse(result["valid"])
        self.assertNotIn("compensate_allele_order\"] =", result["fix"])
        self.assertIn("third way", result["fix"])


class SomalierAlleleOrderUncheckedTest(SimpleTestCase):
    def test_relate_refusing_sites_is_an_answer_not_a_missing_one(self):
        """ Only a somalier before v0.3.5 has no relate --sites, so it can't be the one the settings
            are written for """
        result = allele_order_unchecked("somalier relate failed: Error: unhandled exception: "
                                        "unknown option: --sites [UsageError]")
        self.assertFalse(result["valid"])
        self.assertIn('compensate_allele_order"] = True', result["fix"])

    def test_anything_else_is_a_warning(self):
        result = allele_order_unchecked("tabix failed: no such file")
        self.assertTrue(result["valid"])
        self.assertIn("unchecked", result["warning"])
