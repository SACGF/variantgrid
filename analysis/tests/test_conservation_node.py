from django.test import TestCase, override_settings

from analysis.models.nodes.filters.conservation_node import ConservationNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import create_fake_variant_annotation, create_fake_variants
from annotation.models import VariantAnnotation
from snpdb.models import Variant

GERP = VariantAnnotation.CONSERVATION_SCORES["gerp_pp_rs"]


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestConservationNode(AnalysisSetupMixin, TestCase):
    """ Scaled slider converts a 0-1 fraction into each score's own range """

    def _individual_node(self, **kwargs):
        """ All scores at their minimum (ie no filtering) unless overridden """
        scores = {f: VariantAnnotation.CONSERVATION_SCORES[f]["min"]
                  for f in ConservationNode(analysis=self.analysis).get_individual_field_names()}
        scores.update(kwargs)
        return ConservationNode(analysis=self.analysis, use_individual_sliders=True, **scores)

    def test_scaled_min_zero_does_not_filter(self):
        node = ConservationNode(analysis=self.analysis)
        self.assertFalse(node.modifies_parents())

    def test_scaled_min_above_zero_filters(self):
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.9)
        self.assertTrue(node.modifies_parents())

    def test_scaled_min_uses_each_scores_own_range(self):
        node = ConservationNode(analysis=self.analysis, any_scaled_min=0.9)
        score, _allow_null = node._get_scores_and_allow_null()["gerp_pp_rs"]
        expected = GERP["min"] + 0.9 * (GERP["max"] - GERP["min"])
        self.assertAlmostEqual(score, expected)

    def test_individual_sliders_at_minimum_do_not_filter(self):
        self.assertFalse(self._individual_node().modifies_parents())

    def test_individual_slider_above_minimum_filters(self):
        self.assertTrue(self._individual_node(gerp_pp_rs=4.4).modifies_parents())

    def test_allow_null_skips_scores_that_are_never_null(self):
        """ gerp_pp_rs comes from dbNSFP so intergenic variants are legitimately null """
        node = self._individual_node(gerp_pp_rs=4.4, allow_null=True)
        q_str = str(node._get_node_q())
        self.assertIn("phylop_100_way_vertebrate__isnull", q_str)
        self.assertNotIn("gerp_pp_rs__isnull", q_str)

    def test_filter_returns_expected_rows(self):
        create_fake_variants(self.grch37)
        vav = self.analysis.annotation_version.variant_annotation_version
        variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        v_conserved, v_not_conserved, v_null = variants

        for variant, gerp in [(v_conserved, 5.0), (v_not_conserved, 1.0)]:
            va = create_fake_variant_annotation(variant, vav)
            va.gerp_pp_rs = gerp
            va.save()
        create_fake_variant_annotation(v_null, vav)  # all conservation scores stay null

        node = self._individual_node(gerp_pp_rs=4.4)
        matched = set(Variant.objects.filter(node._get_node_q())
                      .filter(pk__in=[v.pk for v in variants])
                      .values_list("pk", flat=True))
        self.assertEqual(matched, {v_conserved.pk})
