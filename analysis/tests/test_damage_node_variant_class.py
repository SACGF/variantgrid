import itertools

from django.db.models import Q
from django.test import TestCase

from analysis.models.nodes.filters.damage_node import DamageNode, StructuralFilter
from analysis.tests.utils import AnalysisSetupMixin
from annotation.models.damage_enums import PathogenicityImpact
from library.genomics.vcf_enums import VARIANT_CLASS_GROUPS, VariantClass
from snpdb.models import Variant


class VariantClassGroupsTest(TestCase):
    def test_every_variant_class_in_exactly_one_group(self):
        """ The editor only offers what's in a group, so a new VariantClass must be assigned one """
        grouped = list(itertools.chain.from_iterable(VARIANT_CLASS_GROUPS.values()))
        self.assertEqual(len(grouped), len(set(grouped)), "A VariantClass appears in more than one group")
        self.assertEqual(set(grouped), set(VariantClass))


class DamageNodeVariantClassTest(AnalysisSetupMixin, TestCase):
    def test_modifies_parents_with_only_variant_class(self):
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV])
        self.assertTrue(node.modifies_parents())

    def test_empty_variant_class_means_every_type(self):
        node = DamageNode(analysis=self.analysis, variant_class=[])
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_include_restricts_unconditionally(self):
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV])
        self.assertEqual(node._get_node_q(), Q(variantannotation__variant_class__in=[VariantClass.SNV]))

    def test_exclude_negates(self):
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV], variant_class_exclude=True)
        self.assertEqual(node._get_node_q(), ~Q(variantannotation__variant_class__in=[VariantClass.SNV]))

    def test_variant_class_ands_with_scoring_pool(self):
        """ It lands in and_filters, so it restricts rather than joining the damage OR pool """
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV],
                          impact_min=PathogenicityImpact.HIGH)
        q = node._get_node_q()
        self.assertEqual(q.connector, Q.AND)
        self.assertIn(("variantannotation__variant_class__in", [VariantClass.SNV]), q.children)


class DamageNodeStructuralTest(AnalysisSetupMixin, TestCase):
    """ "Just the SVs/CNVs" is a question VEP's variant_class cannot answer - it calls a 1 Mb <DEL>
        and a 1 bp deletion the same class """

    def test_any_adds_nothing(self):
        node = DamageNode(analysis=self.analysis)
        self.assertEqual(StructuralFilter.ANY, node.structural)
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_only_keeps_symbolic(self):
        node = DamageNode(analysis=self.analysis, structural=StructuralFilter.ONLY)
        self.assertTrue(node.modifies_parents())
        self.assertEqual(Variant.get_symbolic_q(), node._get_node_q())

    def test_exclude_negates(self):
        node = DamageNode(analysis=self.analysis, structural=StructuralFilter.EXCLUDE)
        self.assertTrue(node.modifies_parents())
        self.assertEqual(~Variant.get_symbolic_q(), node._get_node_q())

    def test_ands_with_the_class_restriction(self):
        """ Orthogonal filters - both restrict, neither joins the damage OR pool """
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.DELETION],
                          structural=StructuralFilter.ONLY)
        q = node._get_node_q()
        self.assertEqual(Q.AND, q.connector)
        self.assertIn(("variantannotation__variant_class__in", [VariantClass.DELETION]), q.children)
        self.assertIn(("svlen__isnull", False), q.children)
