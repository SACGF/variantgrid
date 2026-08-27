from django.test import TestCase

from annotation.models import VariantAnnotation


class TestAnnotationModels(TestCase):
    def test_variant_annotation_protein_position_to_int(self):
        protein_position_and_expected_int = [
            ("185", 185),
            ("185-187", 185),
            ("?-187", 187),
            ("185-?", 185)
        ]

        for protein_position, expected in protein_position_and_expected_int:
            result = VariantAnnotation.protein_position_to_int(protein_position)
            self.assertEqual(expected, result, f"protein_position_to_int({protein_position})")

    def test_get_short_label_prefers_protein_change(self):
        va = VariantAnnotation(symbol="RUNX1", hgvs_c="NM_001754.5:c.1640C>T",
                               hgvs_p="NP_001745.2:p.Ala547Val")
        self.assertEqual("RUNX1:p.Ala547Val", va.get_short_label())

    def test_get_short_label_falls_back_to_c_hgvs(self):
        va = VariantAnnotation(symbol="RUNX1", hgvs_c="NM_001754.5:c.1640C>T")
        self.assertEqual("RUNX1:c.1640C>T", va.get_short_label())

    def test_get_short_label_without_symbol(self):
        va = VariantAnnotation(hgvs_p="NP_001745.2:p.Ala547Val")
        self.assertEqual("p.Ala547Val", va.get_short_label())

    def test_get_short_label_falls_back_to_g_hgvs(self):
        va = VariantAnnotation(hgvs_g="NC_000021.9:g.34859474G>A")
        self.assertEqual("NC_000021.9:g.34859474G>A", va.get_short_label())
