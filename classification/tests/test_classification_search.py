from django.test import SimpleTestCase

from annotation.models import VariantAnnotation
from classification.signals.classification_search import _variant_location_summary


class VariantLocationSummaryTestCase(SimpleTestCase):
    """ Where an unclassified variant with no canonical c.HGVS lies, shown in its search preview """

    def test_in_gene_without_transcript(self):
        """ RefSeq annotation of MT genes has a symbol but no transcript """
        va = VariantAnnotation(consequence="frameshift_variant", symbol="ND4")
        self.assertEqual(_variant_location_summary(va), "in ND4")

    def test_upstream(self):
        va = VariantAnnotation(consequence="upstream_gene_variant", distance=3709, symbol="CYTB")
        self.assertEqual(_variant_location_summary(va), "3709bp upstream of CYTB")

    def test_intergenic(self):
        va = VariantAnnotation(consequence="intergenic_variant")
        self.assertEqual(_variant_location_summary(va), "intergenic")
