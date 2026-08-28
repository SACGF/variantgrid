from django.test import TestCase

from annotation.models.models_enums import ClinVarOncogenicity, ClinVarPathogenicity
from annotation.vcf_files.clinvar_significance import (
    NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT,
    highest_oncogenicity,
    highest_pathogenicity,
    somatic_tier,
)
from classification.enums import SomaticClinicalSignificance


class SomaticTierTest(TestCase):
    def test_clinvar_tier_strings(self):
        self.assertEqual(somatic_tier("Tier_I_-_Strong"), SomaticClinicalSignificance.TIER_1)
        self.assertEqual(somatic_tier("Tier_II_-_Potential"), SomaticClinicalSignificance.TIER_2)
        self.assertEqual(somatic_tier("Tier_III_-_Unknown"), SomaticClinicalSignificance.TIER_3)
        self.assertEqual(somatic_tier("Tier_IV_-_Benign/Likely_benign"), SomaticClinicalSignificance.TIER_4)

    def test_haplotype_sentinel(self):
        """ The classification belongs to a haplotype, reported through SCIINCL - nothing to derive """
        self.assertIsNone(somatic_tier(NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT))

    def test_no_sci(self):
        self.assertIsNone(somatic_tier(None))


class HighestOncogenicityTest(TestCase):
    def test_case_sensitive_keys_separate_likely_from_plain(self):
        self.assertEqual(highest_oncogenicity("Likely_oncogenic", None), ClinVarOncogenicity.LIKELY_ONCOGENIC)
        self.assertEqual(highest_oncogenicity("Oncogenic", None), ClinVarOncogenicity.ONCOGENIC)
        self.assertEqual(highest_oncogenicity("Likely_benign", None), ClinVarOncogenicity.LIKELY_BENIGN)
        self.assertEqual(highest_oncogenicity("Benign", None), ClinVarOncogenicity.BENIGN)
        self.assertEqual(highest_oncogenicity("Uncertain_significance", None), ClinVarOncogenicity.UNCERTAIN)

    def test_conflicting_reads_oncconf(self):
        conflicting = "Likely_oncogenic(1)|Uncertain_significance(1)"
        self.assertEqual(highest_oncogenicity("Conflicting_classifications_of_oncogenicity", conflicting),
                         ClinVarOncogenicity.LIKELY_ONCOGENIC)

    def test_haplotype_sentinel(self):
        self.assertEqual(highest_oncogenicity(NO_CLASSIFICATION_FOR_THE_SINGLE_VARIANT, None), 0)

    def test_no_onc(self):
        self.assertEqual(highest_oncogenicity(None, None), 0)


class HighestPathogenicityTest(TestCase):
    def test_combined_form(self):
        self.assertEqual(highest_pathogenicity("Benign/Likely_benign", None), ClinVarPathogenicity.LIKELY_BENIGN)

    def test_conflicting_reads_clnsigconf(self):
        conflicting = "Pathogenic(2)|Uncertain_significance(1)"
        self.assertEqual(highest_pathogenicity("Conflicting_interpretations_of_pathogenicity", conflicting),
                         ClinVarPathogenicity.PATHOGENIC)

    def test_conflicting_without_clnsigconf_is_undecidable(self):
        """ None so the importer can count these and die if ClinVar has changed its INFO fields """
        self.assertIsNone(highest_pathogenicity("Conflicting_interpretations_of_pathogenicity", None))

    def test_no_clnsig(self):
        self.assertEqual(highest_pathogenicity(None, None), 0)
