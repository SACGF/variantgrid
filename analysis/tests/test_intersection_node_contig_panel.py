"""
IntersectionNode's Contig panel: whole chromosomes OR cytogenetic bands picked in the form.
"""
from django.test import TestCase

from analysis.forms.forms_nodes import IntersectionNodeForm
from analysis.models.nodes.filters.intersection_node import IntersectionNode
from analysis.tests.utils import AnalysisSetupMixin
from snpdb.models.models_genome import Cytoband
from snpdb.models.models_variant import Variant
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class IntersectionNodeContigPanelTest(AnalysisSetupMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        contig_11 = cls.grch37.chrom_contig_mappings["11"]
        cls.band_11p13 = Cytoband.objects.get(contig=contig_11, name="p13")  # 31000001-36400000
        cls.in_band = slowly_create_test_variant("11", 31000001, "A", "G", cls.grch37)
        cls.outside_band = slowly_create_test_variant("11", 36400001, "A", "G", cls.grch37)
        cls.on_chrom_2 = slowly_create_test_variant("2", 1000, "A", "G", cls.grch37)

    def _save_form(self, contigs=(), cytobands=()) -> IntersectionNode:
        node = IntersectionNode.objects.create(analysis=self.analysis)
        data = {"accordion_panel": IntersectionNode.CONTIG,
                "contigs": [c.pk for c in contigs],
                "cytobands": [c.pk for c in cytobands]}
        form = IntersectionNodeForm(data, instance=node, genome_build=self.grch37)
        self.assertTrue(form.is_valid(), form.errors)
        form.save()
        return IntersectionNode.objects.get(pk=node.pk)

    def _filtered_variants(self, node: IntersectionNode) -> set[Variant]:
        test_variants = [self.in_band, self.outside_band, self.on_chrom_2]
        return set(Variant.objects.filter(node._get_node_q(), pk__in=[v.pk for v in test_variants]))

    def test_band_filters_to_its_coordinates(self):
        node = self._save_form(cytobands=[self.band_11p13])
        self.assertEqual(self._filtered_variants(node), {self.in_band})

    def test_chromosomes_or_bands(self):
        contig_2 = self.grch37.chrom_contig_mappings["2"]
        node = self._save_form(contigs=[contig_2], cytobands=[self.band_11p13])
        self.assertEqual(self._filtered_variants(node), {self.in_band, self.on_chrom_2})

    def test_summary_shows_band_coordinates(self):
        node = self._save_form(cytobands=[self.band_11p13])
        self.assertIn("11p13 = 11:31000001-36400000", node._get_method_summary())
        self.assertEqual(node.get_node_name(), "11p13")
