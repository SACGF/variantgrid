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
