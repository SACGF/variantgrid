""" Node card names and chips say what the node reads (#2599) - the rule is in analysis/AGENTS.md """
from django.test import TestCase

from analysis.enums import PTCNMDFilter, QuadInheritance, ZygosityNodeZygosity
from analysis.models import DamageNode, GeneListNode, QuadNode, SampleNode, ZygosityNode
from analysis.models.nodes.node_display import first_and_more_name_line, truncate_name_line
from analysis.tests.utils import AnalysisSetupMixin
from annotation.enums import PathogenicityImpact
from genes.models import GeneList, SampleGeneList
from snpdb.fake_data import create_fake_quad
from snpdb.models import Sample


class TestNameLineHelpers(TestCase):

    def test_truncate_name_line(self):
        self.assertEqual("x" * 30, truncate_name_line("x" * 30))
        self.assertEqual("x" * 27 + "…", truncate_name_line("x" * 31))

    def test_first_and_more_name_line(self):
        self.assertEqual("", first_and_more_name_line([]))
        self.assertEqual("one", first_and_more_name_line(["one"]))
        line = first_and_more_name_line(["a long gene list name that won't fit", "two", "three"])
        self.assertTrue(line.endswith("… +2 more"), line)
        self.assertLessEqual(len(line), 30)


class TestNodeNames(AnalysisSetupMixin, TestCase):

    def test_quad_node_leads_with_family(self):
        quad = create_fake_quad(self.analysis.user, self.grch37, name="Smith family")
        node = QuadNode(analysis=self.analysis, quad=quad, inheritance=QuadInheritance.RECESSIVE,
                        require_parent_zygosity=True, require_sibling_zygosity=True)
        self.assertEqual(["Smith family", "Recessive"], node.get_node_name().split("\n"))
        node.require_sibling_zygosity = False
        self.assertEqual("Recessive?", node.get_node_name().split("\n")[1])

        quad.name = ""
        self.assertEqual(f"Quad {quad.pk}", node.get_node_name().split("\n")[0])
        node.quad = None
        self.assertEqual("Recessive?", node.get_node_name())

    def test_damage_node_names_headline_filters(self):
        node = DamageNode(analysis=self.analysis, impact_min=PathogenicityImpact.MODERATE,
                          damage_predictions_min=2, ptc_nmd=PTCNMDFilter.ESCAPING)
        lines = node.get_node_name().split("\n")
        self.assertEqual("≥ Moderate impact", lines[0])
        self.assertTrue(lines[1].startswith("2 "), lines[1])
        self.assertEqual("Escapes NMD", lines[2])

    def test_damage_node_score_chips(self):
        node = DamageNode(analysis=self.analysis, popeve_score_max=-5.0, popeve_score_required=True,
                          popeve_score_allow_null=False, revel_score_min=0.7)
        chips = {chip.text: chip.title for chip in node.get_node_chips()}
        self.assertEqual("popEVE score ≤ -5 (required)", chips["popEVE score*"])
        self.assertEqual("REVEL ≥ 0.7 (allows null)", chips["REVEL"])

    def test_gene_list_node_several_lists(self):
        node = GeneListNode.objects.create(analysis=self.analysis)
        for name in ["Cardiomyopathy", "Epilepsy"]:
            gene_list = GeneList.objects.create(name=name, user=self.analysis.user)
            node.genelistnodegenelist_set.create(gene_list=gene_list)
        self.assertEqual("Cardiomyopathy +1 more", node.get_node_name())
        self.assertEqual(["Cardiomyopathy", "Epilepsy"], [chip.title for chip in node.get_node_chips()])

    def test_sample_node_qc_gene_list_chip(self):
        node = SampleNode(analysis=self.analysis, restrict_to_qc_gene_list=True)
        chip = node.get_node_chips()[0]
        self.assertEqual("node-chip-warning", chip.css_class)
        self.assertEqual("Restricted to Sample Gene List, but none specified!", chip.title)

        sample = create_fake_quad(self.analysis.user, self.grch37).proband.sample
        gene_list = GeneList.objects.create(name="QC genes", user=self.analysis.user)
        node.sample_gene_list = SampleGeneList.objects.create(sample=sample, gene_list=gene_list)
        chip = node.get_node_chips()[0]
        self.assertIsNone(chip.css_class)
        self.assertIn("QC genes", chip.title)

    def test_zygosity_node_leads_with_sample(self):
        node = ZygosityNode(analysis=self.analysis, sample=Sample(name="proband"),
                            zygosity=ZygosityNodeZygosity.HET, exclude=True)
        self.assertEqual("proband\nNot HET", node.get_node_name())
        node.zygosity = ZygosityNodeZygosity.MULTIPLE_HIT
        self.assertEqual(["Not Multiple hits in gene"], node.get_node_name().split("\n"))
