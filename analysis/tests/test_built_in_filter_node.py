from django.test import TestCase, override_settings

from analysis.models.nodes.filters.built_in_filter_node import BuiltInFilterNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_annotation import create_fake_variant_annotation, create_fake_variants
from annotation.models.damage_enums import PathogenicityImpact
from snpdb.models import BuiltInFilters, Variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestBuiltInFilterNode(AnalysisSetupMixin, TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        create_fake_variants(cls.grch37)
        vav = cls.analysis.annotation_version.variant_annotation_version
        cls.variants = list(Variant.objects.filter(Variant.get_no_reference_q())[:3])
        cls.va_list = [create_fake_variant_annotation(v, vav) for v in cls.variants]

    def _matched(self, node) -> set:
        return set(Variant.objects.filter(node._get_node_q())
                   .filter(pk__in=[v.pk for v in self.variants])
                   .values_list("pk", flat=True))

    def test_no_filter_selected_does_not_modify_parents(self):
        self.assertFalse(BuiltInFilterNode(analysis=self.analysis).modifies_parents())

    def test_impact_returns_high_and_moderate(self):
        for va, impact in zip(self.va_list, [PathogenicityImpact.HIGH,
                                             PathogenicityImpact.MODERATE,
                                             PathogenicityImpact.LOW]):
            va.impact = impact
            va.save()

        node = BuiltInFilterNode(analysis=self.analysis,
                                 built_in_filter=BuiltInFilters.IMPACT_HIGH_OR_MODERATE)
        v_high, v_moderate, _v_low = self.variants
        self.assertEqual(self._matched(node), {v_high.pk, v_moderate.pk})

    def test_cosmic_returns_variants_with_a_cosmic_id(self):
        va_cosmic = self.va_list[0]
        va_cosmic.cosmic_id = "COSV123"
        va_cosmic.cosmic_count = 5
        va_cosmic.save()

        node = BuiltInFilterNode(analysis=self.analysis, built_in_filter=BuiltInFilters.COSMIC)
        self.assertEqual(self._matched(node), {self.variants[0].pk})

    def test_cosmic_count_min_excludes_rarely_seen_variants(self):
        for va, count in zip(self.va_list, [50, 5, None]):
            va.cosmic_id = "COSV123"
            va.cosmic_count = count
            va.save()

        node = BuiltInFilterNode(analysis=self.analysis, built_in_filter=BuiltInFilters.COSMIC,
                                 cosmic_count_min=10)
        self.assertEqual(self._matched(node), {self.variants[0].pk})

    def test_clinvar_stars_min_adds_review_status_filter(self):
        node = BuiltInFilterNode(analysis=self.analysis, built_in_filter=BuiltInFilters.CLINVAR,
                                 clinvar_stars_min=2)
        self.assertIn("clinvar__review_status__in", str(node._get_node_q()))

    def test_node_name_shows_stars(self):
        node = BuiltInFilterNode(analysis=self.analysis, built_in_filter=BuiltInFilters.CLINVAR,
                                 clinvar_stars_min=3)
        self.assertIn("★★★", node.get_node_name())
