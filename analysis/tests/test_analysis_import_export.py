from django.test import SimpleTestCase

from analysis.analysis_import_export import upgrade_legacy_node_fields
from patients.enums import ZygosityRequirement


class TestLegacyNodeFields(SimpleTestCase):
    """ Exports written before #2659 carry the boolean the ZygosityRequirement choice replaced """

    def test_booleans_become_the_requirement(self):
        data = {"require_parent_zygosity": False, "require_sibling_zygosity": True, "inheritance": "R"}
        upgrade_legacy_node_fields("analysis.QuadNode", data)
        self.assertEqual({"parent_zygosity_requirement": ZygosityRequirement.ALLOW_NO_CALL,
                          "sibling_zygosity_requirement": ZygosityRequirement.CALL_REQUIRED,
                          "inheritance": "R"}, data)

    def test_current_export_is_left_alone(self):
        data = {"parent_zygosity_requirement": ZygosityRequirement.REF_CALL_REQUIRED}
        upgrade_legacy_node_fields("analysis.TrioNode", data)
        self.assertEqual({"parent_zygosity_requirement": ZygosityRequirement.REF_CALL_REQUIRED}, data)
