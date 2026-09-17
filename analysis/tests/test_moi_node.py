from django.test import TestCase, override_settings

from analysis.models import MOINode
from analysis.tests.utils import AnalysisSetupMixin
from patients.models_enums import Zygosity


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestMOINodeZygosities(AnalysisSetupMixin, TestCase):
    """ Each GenCC mode of inheritance maps to the zygosities that can explain it """

    def _zygosities(self, moi, require_zygosity=True):
        node = MOINode(analysis=self.analysis, require_zygosity=require_zygosity)
        return node._get_zygosities(moi)

    def test_dominant_requires_het(self):
        self.assertEqual(self._zygosities("Autosomal dominant"), {Zygosity.HET})

    def test_dominant_without_require_zygosity_allows_hom_alt(self):
        self.assertEqual(self._zygosities("Autosomal dominant", require_zygosity=False),
                         {Zygosity.HET, Zygosity.HOM_ALT})

    def test_x_linked_recessive_requires_hom_alt(self):
        self.assertEqual(self._zygosities("X-linked recessive"), {Zygosity.HOM_ALT})

    def test_somatic_mosaicism_does_not_filter_zygosity(self):
        """ Mosaic variants can be called as anything, even ref """
        self.assertIsNone(self._zygosities("Somatic mosaicism"))

    def test_no_terms_does_not_modify_parents(self):
        node = MOINode.objects.create(analysis=self.analysis)
        self.assertFalse(node.modifies_parents())
