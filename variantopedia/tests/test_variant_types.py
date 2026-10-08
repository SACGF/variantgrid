from django.test import SimpleTestCase, override_settings

from library.genomics import format_bp
from variantopedia.variant_types import get_size_limits, get_unsupported, get_variant_type_rows


@override_settings(VARIANT_SYMBOLIC_ALT_ENABLED=True, VARIANT_SYMBOLIC_ALT_SIZE=50,
                   ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE=1000, ANNOTATION_VEP_SV_MAX_SIZE=None)
class VariantTypesTest(SimpleTestCase):

    def test_del_dup_inv_share_bands_split_at_thresholds(self):
        """ Bands split only where storage or annotation changes (#1358); the inversion's off-by-one is a
            footnote, not a row """
        groups = {g["kind"]: g["bands"] for g in get_variant_type_rows()}
        bands = groups["Deletion, duplication, inversion"]
        self.assertEqual([("Short", "< 50 bp"), ("Small SV", "50 bp to 1 kb"), ("SV", "≥ 1 kb")],
                         [(b["name"], b["range"]) for b in bands])
        self.assertEqual("Sequence (ref/alt)", bands[0]["stored_as"])
        self.assertEqual("Symbolic <DEL>/<DUP>/<INV> + SVLEN", bands[1]["stored_as"])
        self.assertEqual("Standard Short Variant (given to VEP as its sequence)", bands[1]["pipeline"])
        self.assertEqual("Structural Variant", bands[2]["pipeline"])
        self.assertEqual([("Any size", "")], [(b["name"], b["range"]) for b in groups["Insertion, complex substitution"]])

    @override_settings(ANNOTATION_VEP_SV_MAX_SIZE=10_000_000)
    def test_vep_sv_max_size_splits_structural_variants(self):
        """ VEP is never given an SV over its ceiling, so the table must not imply it is annotated as one """
        groups = {g["kind"]: g["bands"] for g in get_variant_type_rows()}
        bands = groups["Deletion, duplication, inversion"]
        self.assertEqual([("SV", "1 kb to 10 Mb"), ("Large SV", "> 10 Mb")],
                         [(b["name"], b["range"]) for b in bands[2:]])
        self.assertIn("not annotated by VEP", bands[3]["pipeline"])
        cnv_bands = groups["Copy number (<CNV>)"]
        self.assertEqual([("SV", "≤ 10 Mb"), ("Large SV", "> 10 Mb")], [(b["name"], b["range"]) for b in cnv_bands])

    @override_settings(VARIANT_SYMBOLIC_ALT_ENABLED=False)
    def test_no_symbolic_storage_is_one_band(self):
        groups = {g["kind"]: g["bands"] for g in get_variant_type_rows()}
        self.assertEqual(["Any size"], [b["name"] for b in groups["Deletion, duplication, inversion"]])

    def test_format_bp(self):
        self.assertEqual(["50 bp", "1 kb", "1.5 kb", "55.12 kb", "10 Mb", "224.23 Mb"],
                         [format_bp(n) for n in (50, 1000, 1500, 55_123, 10_000_000, 224_225_011)])

    @override_settings(LIFTOVER_BCFTOOLS_ENABLED=True, LIFTOVER_BCFTOOLS_MAX_LENGTH=1000, LIFTOVER_BCFTOOLS_SYMBOLIC=False)
    def test_size_limits(self):
        limits = {limit["name"]: limit["text"] for limit in get_size_limits()}
        self.assertIn("Up to 1 kb", limits["Liftover (BCFtools)"])
        self.assertIn("Longer variants are not lifted over", limits["Liftover (BCFtools)"])
        self.assertIn("a duplication up to 5 kb", limits["ClinGen Allele"])
        self.assertIn("<CNV>", limits["g.HGVS"])
        self.assertNotIn("<INS>", limits["g.HGVS"])

    def test_symbolic_kinds_follow_valid_types(self):
        """ <INS> is dropped at import (vcf_clean_alts), so it is listed as unsupported rather than in the table """
        kinds = [g["kind"] for g in get_variant_type_rows()]
        self.assertIn("Copy number (<CNV>)", kinds)
        self.assertNotIn("Insertion (<INS>)", kinds)
        self.assertIn("<INS>", [u["name"] for u in get_unsupported()])

        with override_settings(VARIANT_SYMBOLIC_ALT_VALID_TYPES={"<CNV>", "<DEL>", "<DUP>", "<INV>", "<INS>"}):
            kinds = [g["kind"] for g in get_variant_type_rows()]
            self.assertIn("Copy number <CNV>, insertion <INS>", kinds)
            self.assertNotIn("<INS>", [u["name"] for u in get_unsupported()])
