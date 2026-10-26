"""
The deployment check that keeps settings.SOMALIER["compensate_allele_order"] honest against the
somalier that is actually installed (@see snpdb/variants_to_vcf.py:somalier_alleles_flipped).
"""
from django.test import SimpleTestCase

from variantgrid.deployment_validation.somalier_check import (
    ALLELE_ORDER_SITES_PER_ORDER,
    allele_order_result,
)

NUM_SITES = 2 * ALLELE_ORDER_SITES_PER_ORDER
HALF = ALLELE_ORDER_SITES_PER_ORDER


class SomalierAlleleOrderResultTest(SimpleTestCase):
    """ Every site the check writes is hom-alt against a no-call, so how many come back as opposite
        homozygotes - uncompensated, then compensated if that was wrong - says which way to set it """

    def test_fixed_somalier_must_not_be_compensated(self):
        """ Compensating a somalier that reads REF/ALT itself breaks ancestry, not relate, so only the
            uncompensated probe can tell """
        self.assertTrue(allele_order_result("GRCh38", NUM_SITES, NUM_SITES, None, False)["valid"])
        result = allele_order_result("GRCh38", NUM_SITES, NUM_SITES, None, True)
        self.assertFalse(result["valid"])
        self.assertIn('compensate_allele_order"] = False', result["fix"])
        self.assertIn("ancestry", result["fix"])
        self.assertIn("somalier_existing_vcfs --clear", result["fix"])

    def test_unfixed_somalier_must_be_compensated(self):
        """ v0.3.5 (half the sites wrong uncompensated) and older (no relate --sites at all) """
        for uncompensated_ibs0 in (HALF, None):
            with self.subTest(uncompensated_ibs0=uncompensated_ibs0):
                self.assertTrue(allele_order_result("GRCh38", NUM_SITES, uncompensated_ibs0, NUM_SITES, True)["valid"])
                result = allele_order_result("GRCh38", NUM_SITES, uncompensated_ibs0, NUM_SITES, False)
                self.assertFalse(result["valid"])
                self.assertIn('compensate_allele_order"] = True', result["fix"])

    def test_neither_way_right_does_not_claim_to_know_why(self):
        for compensating in (True, False):
            with self.subTest(compensating=compensating):
                result = allele_order_result("GRCh38", NUM_SITES, HALF, 0, compensating)
                self.assertFalse(result["valid"])
                self.assertNotIn('compensate_allele_order"] =', result["fix"])
                self.assertIn("third way", result["fix"])
