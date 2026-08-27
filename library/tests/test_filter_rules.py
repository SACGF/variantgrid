from django.test import SimpleTestCase

from library.django_utils.filter_rules import filter_rules_to_q


def _q(op, field="variantannotation__dbsnp_rs_id", data="", group_op="AND", **kwargs):
    return filter_rules_to_q({"groupOp": group_op, "rules": [{"op": op, "field": field, "data": data}]}, **kwargs)


class FilterRulesToQTest(SimpleTestCase):
    """ The rule -> Django lookup translation shared by the grids' filter builder and FilterNode """

    def test_is_null_is_not_an_inverted_q(self):
        """ ~Q(field__isnull=True) generated a full table scan on variant (>100M rows) - the negation
            goes in as the lookup value instead """
        self.assertNotIn("NOT", str(_q("nu")))
        self.assertIn("('variantannotation__dbsnp_rs_id__isnull', True)", str(_q("nu")))
        self.assertIn("('variantannotation__dbsnp_rs_id__isnull', False)", str(_q("nn")))

    def test_in_splits_on_commas(self):
        q = _q("in", field="locus__contig__name", data="1,2,X")
        self.assertIn("['1', '2', 'X']", str(q))

    def test_negated_operations_invert(self):
        self.assertIn("NOT", str(_q("ne", data="rs123")))
        self.assertNotIn("NOT", str(_q("eq", data="rs123")))

    def test_ignore_case_by_default(self):
        self.assertIn("icontains", str(_q("cn", data="foo")))
        self.assertIn("contains", str(_q("cn", data="foo", ignore_case=False)))
        self.assertNotIn("icontains", str(_q("cn", data="foo", ignore_case=False)))

    def test_group_op(self):
        rules = [{"op": "gt", "field": "locus__position", "data": "10"},
                 {"op": "lt", "field": "locus__position", "data": "20"}]
        self.assertIn("OR", str(filter_rules_to_q({"groupOp": "OR", "rules": rules})))
        self.assertIn("AND", str(filter_rules_to_q({"groupOp": "AND", "rules": rules})))

    def test_no_rules(self):
        self.assertIsNone(filter_rules_to_q(None))
        self.assertIsNone(filter_rules_to_q({"groupOp": "AND", "rules": []}))
