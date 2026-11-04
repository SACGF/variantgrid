"""The Gene Symbol Matching page: the rule it reports for each kind of name, and the typed-name form."""
from collections import defaultdict
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from genes.gene_symbol_matching import GeneNameRule, explain_gene_names
from genes.tests.gene_symbol_test_utils import create_gene_symbol_alias_fixture
from variantgrid.perm_path import get_visible_url_names

EXPECTED_RULES = {
    "SEPTIN2": GeneNameRule.CURRENT_SYMBOL,
    "AIP": GeneNameRule.CURRENT_SYMBOL,  # another gene's alias, but an approved symbol itself
    "sept2": GeneNameRule.PREVIOUS_SYMBOL,  # case-insensitive, and the rename beats SEPTIN6's alias
    "ACPP": GeneNameRule.PREVIOUS_SYMBOL,  # has its own GeneSymbol row, still a rename
    "SEPTX": GeneNameRule.AMBIGUOUS,
    "RETIRED_ONLY": GeneNameRule.ALIAS_WITHOUT_HGNC,
    "NOT_A_GENE_SYMBOL": GeneNameRule.UNKNOWN,
}


class TestGeneSymbolMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_gene_symbol_alias_fixture()
        cls.user = User.objects.create(username="gene_symbol_matching_user")

    def test_rule_for_each_kind_of_name(self):
        with self.assertNumQueries(4):  # resolve_names' three plus the aliases followed
            explanations = explain_gene_names(EXPECTED_RULES)
        for name, rule in EXPECTED_RULES.items():
            self.assertEqual(rule, explanations[name].rule, name)
        self.assertEqual("SEPTIN2", explanations["sept2"].gene_symbol_id)
        self.assertEqual("SEPT2", explanations["sept2"].alias.alias)
        self.assertEqual(["SEPTIN2", "SEPTIN6"], explanations["SEPTX"].candidates)
        self.assertFalse(explanations["NOT_A_GENE_SYMBOL"].matched)

    def test_page_explains_typed_name(self):
        self.client.force_login(self.user)
        response = self.client.get(reverse("gene_symbol_matching"), {"name": " sept2 "})
        self.assertEqual(200, response.status_code)
        self.assertEqual("sept2", response.context["typed_name"])
        self.assertEqual(GeneNameRule.PREVIOUS_SYMBOL, response.context["typed"].rule)

    @override_settings(VARIANT_GENE_LEVEL_ENABLED=False)
    def test_page_leaves_out_features_the_deployment_hides(self):
        """ Shariant has no gene lists or gene level variants, so the page mustn't send users looking for them """
        url_name_visible = defaultdict(lambda: False, get_visible_url_names(), gene_lists=False)
        self.client.force_login(self.user)
        with patch("snpdb.processors.get_visible_url_names", return_value=url_name_visible):
            response = self.client.get(reverse("gene_symbol_matching"))
        for text in ["a gene list,", "Gene lists,", "gene list grid", "PanelApp", "fusions"]:
            self.assertNotContains(response, text)
