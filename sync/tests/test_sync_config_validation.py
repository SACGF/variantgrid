from django.test import SimpleTestCase

from sync.shariant.query_json_filter import QueryJsonFilter
from sync.shariant.variant_grid_download import LAB_GROUP_NAME_PATTERN, _group_names_param


class QueryJsonFilterValidationTest(SimpleTestCase):

    def test_evidence_keys_accepted(self):
        json_filter = QueryJsonFilter.classification_value_filter()
        json_filter.convert_to_q({"acmg:pvs1": "PM", "1000_genomes_af": {"lt": 0.1}, "allele_origin": ["germline"]})

    def test_extra_lookups_rejected(self):
        json_filter = QueryJsonFilter.classification_value_filter()
        for blob in [{"allele_origin__isnull": True}, {"allele_origin": {"value__isnull": True}}]:
            with self.assertRaises(ValueError, msg=str(blob)):
                json_filter.convert_to_q(blob)


class DownloadGroupNamesTest(SimpleTestCase):

    def test_group_names_param(self):
        config = {"exclude_labs": ["org_a/lab-1", "org_b/lab_2"]}
        self.assertEqual(_group_names_param(config, "exclude_labs", LAB_GROUP_NAME_PATTERN), "org_a/lab-1,org_b/lab_2")
        self.assertIsNone(_group_names_param({}, "exclude_labs", LAB_GROUP_NAME_PATTERN))

    def test_malformed_group_name_rejected(self):
        for group_name in ["org_a", "org_a/lab,org_b/lab", "org_a/lab/extra"]:
            with self.assertRaises(ValueError, msg=group_name):
                _group_names_param({"exclude_labs": [group_name]}, "exclude_labs", LAB_GROUP_NAME_PATTERN)
