from django.test import TestCase, override_settings

from analysis.models import GroupOperation
from analysis.models.nodes.filters.tissue_node import TissueNode
from analysis.tests.utils import AnalysisSetupMixin


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTissueNodeText(AnalysisSetupMixin, TestCase):
    """ UniProtKB panel matches words against the tissue specificity text """

    def _node(self, text_tissue, group_operation=GroupOperation.ALL):
        return TissueNode(analysis=self.analysis, accordion_panel=TissueNode.UNIPROTKB,
                          text_tissue=text_tissue, group_operation=group_operation)

    def test_no_text_does_not_modify_parents(self):
        self.assertFalse(self._node(None).modifies_parents())

    def test_text_modifies_parents(self):
        self.assertTrue(self._node("kidney").modifies_parents())

    def test_hpa_panel_needs_a_tissue_sample(self):
        node = TissueNode(analysis=self.analysis, accordion_panel=TissueNode.HUMAN_PROTEIN_ATLAS,
                          text_tissue="kidney")
        self.assertFalse(node.modifies_parents())

    def test_all_requires_every_word(self):
        q_str = str(self._node("kidney liver", group_operation=GroupOperation.ALL)._get_node_q())
        self.assertIn("AND", q_str)
        self.assertIn("kidney", q_str)
        self.assertIn("liver", q_str)

    def test_any_matches_either_word(self):
        q_str = str(self._node("kidney liver", group_operation=GroupOperation.ANY)._get_node_q())
        self.assertIn("OR", q_str)
        self.assertIn("kidney", q_str)
        self.assertIn("liver", q_str)
