from django.test import SimpleTestCase, override_settings

from variantopedia.variant_types import VARIANT_TYPE_COLUMNS, get_variant_type_rows


class VariantTypesTest(SimpleTestCase):

    @override_settings(VARIANT_SYMBOLIC_ALT_SIZE=50, ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE=1000)
    def test_deletion_bands_follow_thresholds(self):
        """ Bands split where a rule changes (#1358) """
        stored_as = VARIANT_TYPE_COLUMNS.index("Stored as")
        pipeline = VARIANT_TYPE_COLUMNS.index("Annotation pipeline")
        deletions = {r["size"]: r["cells"] for r in get_variant_type_rows() if r["kind"] == "Deletion"}
        self.assertEqual("Sequence (ref/alt)", deletions["1-49bp"][stored_as])
        self.assertEqual(("Symbolic <DEL> + SVLEN", "Standard Short Variant"),
                         (deletions["50-999bp"][stored_as], deletions["50-999bp"][pipeline]))
        self.assertEqual("Structural Variant", deletions["1,000-10,000bp"][pipeline])
