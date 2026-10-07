from django.test import SimpleTestCase, override_settings

from variantopedia.variant_types import format_bp, get_size_limits, get_variant_type_rows


@override_settings(VARIANT_SYMBOLIC_ALT_ENABLED=True, VARIANT_SYMBOLIC_ALT_SIZE=50,
                   ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE=1000)
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

    @override_settings(VARIANT_SYMBOLIC_ALT_ENABLED=False)
    def test_no_symbolic_storage_is_one_band(self):
        groups = {g["kind"]: g["bands"] for g in get_variant_type_rows()}
        self.assertEqual(["Any size"], [b["name"] for b in groups["Deletion, duplication, inversion"]])

    def test_format_bp(self):
        self.assertEqual(["50 bp", "1 kb", "1.5 kb", "10 Mb"], [format_bp(n) for n in (50, 1000, 1500, 10_000_000)])

    @override_settings(LIFTOVER_BCFTOOLS_ENABLED=True, LIFTOVER_BCFTOOLS_MAX_LENGTH=1000, LIFTOVER_BCFTOOLS_SYMBOLIC=False)
    def test_size_limits(self):
        limits = {limit["name"]: limit["text"] for limit in get_size_limits()}
        self.assertIn("Up to 1 kb", limits["Liftover (BCFtools)"])
        self.assertIn("Longer variants are not lifted over", limits["Liftover (BCFtools)"])
        self.assertIn("a duplication up to 5 kb", limits["ClinGen Allele"])
        self.assertIn("<CNV>, <INS>", limits["g.HGVS"])
