from django.test import SimpleTestCase

from library.utils.export_utils import ExportRow, ExportTweak, export_column


class CategorisedRow(ExportRow):

    @export_column()
    def always(self):
        return None

    @export_column(categories={"transient": True})
    def transient(self):
        return None

    @export_column(categories={"count": 0})
    def count_zero(self):
        return None

    @export_column(categories={"system": {"generic", "varseq"}})
    def generic_or_varseq(self):
        return None


class ExportCategoriesTest(SimpleTestCase):

    @staticmethod
    def _columns(categories) -> list[str]:
        export_methods = CategorisedRow.get_export_methods(export_tweak=ExportTweak(categories=categories))
        return [export_method.__name__ for export_method in export_methods]

    def test_unprovided_category_does_not_restrict(self):
        self.assertEqual(self._columns({"other": None}),
                         ["always", "transient", "count_zero", "generic_or_varseq"])

    def test_falsy_provided_values_are_honoured(self):
        self.assertNotIn("transient", self._columns({"transient": None}))
        self.assertNotIn("transient", self._columns({"transient": False}))
        self.assertIn("count_zero", self._columns({"count": 0}))
        self.assertNotIn("count_zero", self._columns({"count": 1}))

    def test_any_tagged_value_matches(self):
        self.assertIn("generic_or_varseq", self._columns({"system": "varseq"}))
        self.assertIn("generic_or_varseq", self._columns({"system": ("core", "generic")}))
        self.assertNotIn("generic_or_varseq", self._columns({"system": "emedgene"}))
