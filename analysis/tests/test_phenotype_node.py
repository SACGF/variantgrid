from django.test import TestCase, override_settings

from analysis.models import PhenotypeNode
from analysis.tests.utils import AnalysisSetupMixin

TEXT_COLUMNS = [
    "variantannotation__gene__summary",
    "variantannotation__gene__geneannotation__omim_terms",
    "variantannotation__transcript_version__gene_version__hgnc__uniprot__function",
    "variantannotation__gene__geneannotation__hpo_terms",
]


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPhenotypeNodeText(AnalysisSetupMixin, TestCase):
    """ Free text is matched against gene descriptions as well as the ontology terms """

    def _node(self, text_phenotype=None) -> PhenotypeNode:
        return PhenotypeNode.objects.create(analysis=self.analysis, text_phenotype=text_phenotype)

    def test_no_terms_or_text_does_not_modify_parents(self):
        self.assertFalse(self._node().modifies_parents())

    def test_text_modifies_parents(self):
        self.assertTrue(self._node("epilepsy").modifies_parents())

    def test_text_searches_every_description_column(self):
        q_str = str(self._node("epilepsy")._get_node_q())
        for column in TEXT_COLUMNS:
            self.assertIn(f"{column}__icontains", q_str)

    def test_each_word_is_searched_separately(self):
        q_str = str(self._node("epilepsy ataxia")._get_node_q())
        self.assertIn("epilepsy", q_str)
        self.assertIn("ataxia", q_str)

    def test_node_name_shows_the_text(self):
        self.assertIn("epilepsy", self._node("epilepsy").get_node_name())
