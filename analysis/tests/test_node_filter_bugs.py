"""
Tests for analysis node filters where a count of 0 is a real user choice rather
than "unset" - these all used to be lost to truthiness checks.
"""

from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode, Analysis
from analysis.models.enums import GroupOperation
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.tests.utils import AnalysisSetupMixin

# ---------------------------------------------------------------------------
# AbstractZygosityCountNode - max_count=0
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestZygosityCountZero(AnalysisSetupMixin, TestCase):
    """ max_count=0 means "seen in no samples" - a real filter, not an unset one. """

    def _node_with_max_zero(self):
        node = AllVariantsNode(analysis=self.analysis)
        node.max_het_or_hom_count = 0
        return node

    def test_max_count_zero_produces_lte_filter(self):
        """max_count=0 must generate a col__lte=0 filter."""
        node = self._node_with_max_zero()
        arg_q_dict = node.get_zygosity_count_arg_q_dict()

        all_q_strs = [
            str(q)
            for q_dict in arg_q_dict.values()
            for q in q_dict.values()
        ]
        self.assertTrue(
            any("__lte" in s for s in all_q_strs),
            f"max_count=0 produced no __lte filter. Got Q strings: {all_q_strs}",
        )

    def test_max_count_zero_differs_from_no_max(self):
        """max_count=0 and max_count=None must produce different Q dicts."""
        node_zero = self._node_with_max_zero()
        node_none = AllVariantsNode(analysis=self.analysis)
        # max_het_or_hom_count defaults to None (nullable field)

        q_zero = str(node_zero.get_zygosity_count_arg_q_dict())
        q_none = str(node_none.get_zygosity_count_arg_q_dict())

        self.assertNotEqual(
            q_zero,
            q_none,
            "max_count=0 and max_count=None produce identical Q dicts",
        )

    def test_description_shows_max_zero(self):
        """_get_zygosity_count_description() must include '0' when max_count=0."""
        node = self._node_with_max_zero()
        desc = node._get_zygosity_count_description()
        self.assertIn("0", desc, f"max_count=0 disappeared from description. Got: {desc!r}")


# ---------------------------------------------------------------------------
# PopulationNode.modifies_parents() - gnomad_hom_alt_max=0
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeModifiesParents(AnalysisSetupMixin, TestCase):
    """
    gnomad_hom_alt_max=0 (remove variants with ANY gnomAD hom-alt carrier) is a real
    filter, not "unset" - modifies_parents() has to report it or _get_node_q()'s Q is
    never applied to the queryset.
    """

    def _node(self, **kwargs):
        """Return an unsaved PopulationNode with 100% (no population AF filter) by default."""
        defaults = dict(percent=PopulationNode.EVERYTHING, show_gnomad_filtered=True)
        defaults.update(kwargs)
        return PopulationNode(analysis=self.analysis, **defaults)

    def test_gnomad_hom_alt_max_zero_modifies_parents(self):
        node = self._node(gnomad_hom_alt_max=0)
        self.assertTrue(node.modifies_parents())

    def test_gnomad_hom_alt_max_zero_q_filter_is_not_none(self):
        node = self._node(gnomad_hom_alt_max=0)
        self.assertIsNotNone(node._get_node_q(),
                             "_get_node_q() returned None for gnomad_hom_alt_max=0")

    def test_gnomad_hom_alt_max_none_does_not_modify_alone(self):
        node = self._node(gnomad_hom_alt_max=None)
        self.assertFalse(node.modifies_parents())

    def test_gnomad_hom_alt_max_positive_modifies_parents(self):
        node = self._node(gnomad_hom_alt_max=5)
        self.assertTrue(node.modifies_parents())


# ---------------------------------------------------------------------------
# PopulationNode group operation Q structure
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeGroupOperationQ(AnalysisSetupMixin, TestCase):
    """
    PopulationNode uses *inverted* operators for ANY/ALL (see comment in source):

        OPERATIONS = {
            GroupOperation.ALL: operator.or_,   # inverted
            GroupOperation.ANY: operator.and_,  # inverted
        }

    Semantics:
      ANY (strict)  — remove if above cutoff in ANY database → AND all per-db filters.
      ALL (lenient) — remove only if above cutoff in ALL databases → OR per-db filters.

    These tests pin the observable Q structure so the inversion isn't "corrected" away.
    """

    def _pop_node_q(self, group_operation):
        """Create a saved node using gnomad_af + af_1kg + af_uk10k at 1% and return _get_node_q().
        Must be saved so that self.populationnodegnomadpopulation_set.all() can run.
        """
        node = PopulationNode.objects.create(
            analysis=self.analysis,
            percent=1.0,
            group_operation=group_operation,
            gnomad_af=True,
            gnomad_popmax_af=False,
            af_1kg=True,
            af_uk10k=True,
            topmed_af=False,
        )
        return node._get_node_q()

    def test_any_group_operation_combines_with_and(self):
        """GroupOperation.ANY → AND of per-db filters (all must pass → strict)."""
        q = self._pop_node_q(GroupOperation.ANY)
        self.assertIsNotNone(q)
        q_str = str(q)
        self.assertTrue(
            q_str.startswith("(AND:"),
            f"ANY mode should combine filters with AND, got: {q_str!r}",
        )

    def test_all_group_operation_combines_with_or(self):
        """GroupOperation.ALL → OR of per-db filters (any can pass → lenient)."""
        q = self._pop_node_q(GroupOperation.ALL)
        self.assertIsNotNone(q)
        q_str = str(q)
        self.assertTrue(
            q_str.startswith("(OR:"),
            f"ALL mode should combine filters with OR, got: {q_str!r}",
        )
