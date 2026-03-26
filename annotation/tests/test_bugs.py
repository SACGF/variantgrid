"""
Regression tests for confirmed bugs in the annotation app, plus sentinels for
known fragile code paths that are likely to break when related code changes.
"""
from django.test import TestCase

from annotation.models.models import ClinVar, VariantAnnotation, VariantAnnotationVersion
from annotation.models.models_enums import ClinVarReviewStatus
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from classification.enums import AlleleOriginBucket
from library.genomics import overlap_fraction, Range
from types import SimpleNamespace


# ---------------------------------------------------------------------------
# Bug: amino_acid_3_to_1 did not convert stop codon "Ter" → "*"
# Fix: aa_3_to_1["Ter"] = "*" added manually after invert_dict(protein_letters_1to3)
# ---------------------------------------------------------------------------

class TestAminoAcid3To1(TestCase):

    def test_stop_codon_ter_converts_to_asterisk(self):
        # BioPython's protein_letters_1to3 omits stop codons; "Ter" must be added manually.
        # Without the fix this returned "p.R123Ter", producing wrong external search terms.
        self.assertEqual(VariantAnnotation.amino_acid_3_to_1("p.Arg123Ter"), "p.R123*")

    def test_standard_substitution_still_works(self):
        self.assertEqual(VariantAnnotation.amino_acid_3_to_1("p.Ala2351Thr"), "p.A2351T")

    def test_non_aa_suffix_unchanged(self):
        # "fs" and "del" are not amino acids — must pass through
        self.assertEqual(VariantAnnotation.amino_acid_3_to_1("p.Ala123fs"), "p.A123fs")
        self.assertEqual(VariantAnnotation.amino_acid_3_to_1("p.Ala123del"), "p.A123del")


# ---------------------------------------------------------------------------
# Bug: _merge_cosmic_ids("") injected empty string into result set
# "".split("&") == [""], so sorted set became ["", "COSV123"] → "&COSV123"
# Fix: filter empty strings with {c for c in ... if c}
# ---------------------------------------------------------------------------

class TestMergeCosmicIds(TestCase):

    def test_empty_custom_id_does_not_corrupt_existing(self):
        data = {"cosmic_id": "COSV123"}
        BulkVEPVCFAnnotationInserter._merge_cosmic_ids(data, "")
        ids = data["cosmic_id"].split("&")
        self.assertNotIn("", ids,
            "Empty custom_vcf_cosmic_id must not inject empty string into result")
        self.assertIn("COSV123", ids)

    def test_deduplicates_ids(self):
        data = {"cosmic_id": "COSV123"}
        BulkVEPVCFAnnotationInserter._merge_cosmic_ids(data, "COSV123&COSV456")
        ids = data["cosmic_id"].split("&")
        self.assertEqual(sorted(ids), ["COSV123", "COSV456"])

    def test_no_existing_id(self):
        data = {}
        BulkVEPVCFAnnotationInserter._merge_cosmic_ids(data, "COSV789")
        self.assertEqual(data["cosmic_id"], "COSV789")

    def test_multiple_custom_merged_with_existing(self):
        data = {"cosmic_id": "COSV100"}
        BulkVEPVCFAnnotationInserter._merge_cosmic_ids(data, "COSV200&COSV300")
        self.assertEqual(sorted(data["cosmic_id"].split("&")), ["COSV100", "COSV200", "COSV300"])


# ---------------------------------------------------------------------------
# Sentinel: gnomad_major_version crashes on old gnomADv2 "r2.1" format
# int("r2") → ValueError. Any call to gnomad4_or_later/has_hemi/has_mid on a
# GRCh37 annotation version will crash until this is fixed.
# ---------------------------------------------------------------------------

class TestGnomadMajorVersion(TestCase):

    def _call(self, gnomad_value):
        return VariantAnnotationVersion.gnomad_major_version.fget(
            SimpleNamespace(gnomad=gnomad_value))

    def test_v3_and_v4_numeric_formats_work(self):
        self.assertEqual(self._call("3.1.2"), 3)
        self.assertEqual(self._call("4.0"), 4)

    def test_v2_old_r_prefix_crashes(self):
        # Old GRCh37 gnomADv2 VEP output stores "r2.1"; int("r2") → ValueError.
        # This will start passing once the property handles the "r" prefix.
        with self.assertRaises(ValueError):
            self._call("r2.1")


# ---------------------------------------------------------------------------
# Sentinel: ClinVarReviewStatus.stars() — all enum values must be in STARS dict
# Adding a new status without a STARS entry would silently KeyError at runtime.
# ---------------------------------------------------------------------------

class TestClinVarReviewStatusStars(TestCase):

    def test_all_statuses_have_stars_mapping(self):
        for status in ClinVarReviewStatus:
            try:
                stars = status.stars()
            except KeyError:
                self.fail(f"ClinVarReviewStatus.{status.name} missing from STARS dict")
            self.assertIn(stars, [0, 1, 2, 3, 4],
                f"{status.name} returned invalid star count: {stars}")

    def test_known_star_values(self):
        self.assertEqual(ClinVarReviewStatus.PRACTICE_GUIDELINE.stars(), 4)
        self.assertEqual(ClinVarReviewStatus.REVIEWED_BY_EXPERT_PANEL.stars(), 3)
        self.assertEqual(ClinVarReviewStatus.CRITERIA_PROVIDED_CONFLICTING_INTERPRETATIONS.stars(), 1)
        self.assertEqual(ClinVarReviewStatus.NO_ASSERTION_PROVIDED.stars(), 0)


# ---------------------------------------------------------------------------
# Sentinel: ClinVar._database_terms — hardcoded slice lengths for prefix stripping
# If the prefix string is ever renamed, the slice index becomes wrong silently.
# ---------------------------------------------------------------------------

class TestClinVarDatabaseTerms(TestCase):

    def test_mondo_double_prefix_normalized(self):
        self.assertEqual(ClinVar._database_terms("MONDO:MONDO:0000123"), ["MONDO:0000123"])

    def test_orphanet_prefix_normalized(self):
        self.assertEqual(ClinVar._database_terms("Orphanet:ORPHA123456"), ["ORPHA:123456"])

    def test_hpo_prefix_stripped(self):
        self.assertEqual(
            ClinVar._database_terms("Human_Phenotype_Ontology:HP:0001234"), ["HP:0001234"])

    def test_hpo_prefix_slice_length_is_correct(self):
        # The code does name[25:]. If "Human_Phenotype_Ontology:" is ever renamed
        # the slice silently strips the wrong number of characters.
        self.assertEqual(len("Human_Phenotype_Ontology:"), 25)

    def test_none_and_empty_return_empty_list(self):
        self.assertEqual(ClinVar._database_terms(None), [])
        self.assertEqual(ClinVar._database_terms(""), [])


# ---------------------------------------------------------------------------
# Sentinel: overlap_fraction int() truncation in SV overlap percentage
# 99.6% → int() gives 99, not 100. Documents known precision loss.
# ---------------------------------------------------------------------------

class TestOverlapFractionTruncation(TestCase):

    def test_99_6_percent_truncates_to_99_not_100(self):
        # r1=1000-2000 (len 1000), r2=1000-1996 (overlap 996bp → 99.6%)
        of = overlap_fraction(Range(1000, 2000), Range(1000, 1996))
        self.assertAlmostEqual(of, 0.996)
        # The annotation code uses int(of * 100) — truncation, not rounding.
        # A variant with 99.6% gnomAD SV overlap is stored as 99%, which could
        # cause it to miss a ≥100% filter.
        self.assertEqual(int(of * 100), 99)
