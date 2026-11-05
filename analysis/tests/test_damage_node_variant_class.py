import itertools

from django.test import SimpleTestCase, TestCase

from analysis.enums import StructuralFilter
from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import AnalysisSetupMixin, annotate_variants, node_keeps
from annotation.enums import PathogenicityImpact
from library.genomics.vcf_enums import VARIANT_CLASS_GROUPS, VariantClass
from snpdb.models import Variant
from snpdb.models.models_variant import VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)


class _AnnotatedVariantsMixin(AnalysisSetupMixin):
    """ Named variants, each with a VariantAnnotation row in the analysis's partition, to filter with a node """

    @classmethod
    def _annotate(cls, variants_and_columns: dict[str, tuple[Variant, dict]]):
        cls.variants = annotate_variants(cls.analysis.annotation_version.variant_annotation_version,
                                         variants_and_columns)

    def _filtered(self, node: DamageNode) -> set[str]:
        return node_keeps(node, self.variants)


class VariantClassGroupsTest(TestCase):
    def test_every_variant_class_in_exactly_one_group(self):
        """ The editor only offers what's in a group, so a new VariantClass must be assigned one """
        grouped = list(itertools.chain.from_iterable(VARIANT_CLASS_GROUPS.values()))
        self.assertEqual(len(grouped), len(set(grouped)), "A VariantClass appears in more than one group")
        self.assertEqual(set(grouped), set(VariantClass))


class DamageNodeVariantClassTest(_AnnotatedVariantsMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        snv_high = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        snv_low = slowly_create_test_variant("1", 1100, "A", "T", cls.grch37)
        insertion_high = slowly_create_test_variant("1", 1200, "A", "AT", cls.grch37)
        null_class = slowly_create_test_variant("1", 1300, "A", "T", cls.grch37)
        cls._annotate({
            "snv_high": (snv_high, {"variant_class": VariantClass.SNV, "impact": PathogenicityImpact.HIGH}),
            "snv_low": (snv_low, {"variant_class": VariantClass.SNV, "impact": PathogenicityImpact.LOW}),
            "insertion_high": (insertion_high, {"variant_class": VariantClass.INSERTION,
                                                "impact": PathogenicityImpact.HIGH}),
            "null_class": (null_class, {"variant_class": None}),
        })

    def test_modifies_parents_with_only_variant_class(self):
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV])
        self.assertTrue(node.modifies_parents())

    def test_empty_variant_class_means_every_type(self):
        node = DamageNode(analysis=self.analysis, variant_class=[])
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_include_restricts_unconditionally(self):
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV])
        self.assertEqual({"snv_high", "snv_low"}, self._filtered(node))

    def test_exclude_negates(self):
        """ A variant VEP gave no class isn't known to be an SNV, so "not SNV" keeps it """
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV], variant_class_exclude=True)
        self.assertEqual({"insertion_high", "null_class"}, self._filtered(node))

    def test_variant_class_ands_with_scoring_pool(self):
        """ It restricts rather than joining the damage OR pool """
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.SNV],
                          impact_min=PathogenicityImpact.HIGH)
        self.assertEqual({"snv_high"}, self._filtered(node))


class DamageNodeTypeChipsTest(SimpleTestCase):
    """ One card chip per variant class group the node restricts to """

    @staticmethod
    def _chip_texts(**kwargs) -> list[str]:
        return [chip.text for chip in DamageNode(**kwargs).get_node_chips()]

    def test_no_restriction_no_chips(self):
        self.assertEqual([], self._chip_texts(variant_class=[]))

    def test_whole_groups_named_for_the_group(self):
        indel = [vc.value for vc in VARIANT_CLASS_GROUPS["Indel"]]
        self.assertEqual(["SNV"], self._chip_texts(variant_class=[VariantClass.SNV.value]))
        self.assertEqual(["SNV", "Indel"], self._chip_texts(variant_class=[*indel, VariantClass.SNV.value]))

    def test_partial_group_lists_its_classes(self):
        chips = DamageNode(variant_class=[VariantClass.DELETION.value, VariantClass.INSERTION.value]).get_node_chips()
        self.assertEqual(["insertion, deletion"], [chip.text for chip in chips])
        self.assertIn("2 of 6 Indel", chips[0].title)

    def test_exclude(self):
        texts = self._chip_texts(variant_class=[VariantClass.SNV.value, VariantClass.DELETION.value],
                                 variant_class_exclude=True)
        self.assertEqual(["not SNV", "not deletion"], texts)


class DamageNodeStructuralTest(_AnnotatedVariantsMixin, TestCase):
    """ "Just the SVs/CNVs" is a question VEP's variant_class cannot answer - it calls a 1 Mb <DEL>
        and a 1 bp deletion the same class """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        sv_deletion = VariantCoordinate(chrom="1", position=10_000, ref="A", alt="<DEL>",
                                        end=1_010_000, svlen=-1_000_000)
        sv_duplication = VariantCoordinate(chrom="1", position=2_000_000, ref="A", alt="<DUP>",
                                           end=2_100_000, svlen=100_000)
        cls._annotate({
            "sv_deletion": (slowly_create_test_variant_from_coordinate(sv_deletion, cls.grch37),
                            {"variant_class": VariantClass.DELETION}),
            "small_deletion": (slowly_create_test_variant("1", 3_000_000, "AT", "A", cls.grch37),
                               {"variant_class": VariantClass.DELETION}),
            "sv_duplication": (slowly_create_test_variant_from_coordinate(sv_duplication, cls.grch37),
                               {"variant_class": VariantClass.DUPLICATION}),
        })

    def test_any_adds_nothing(self):
        node = DamageNode(analysis=self.analysis)
        self.assertEqual(StructuralFilter.ANY, node.structural)
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_only_keeps_symbolic(self):
        node = DamageNode(analysis=self.analysis, structural=StructuralFilter.ONLY)
        self.assertTrue(node.modifies_parents())
        self.assertEqual({"sv_deletion", "sv_duplication"}, self._filtered(node))

    def test_exclude_negates(self):
        node = DamageNode(analysis=self.analysis, structural=StructuralFilter.EXCLUDE)
        self.assertTrue(node.modifies_parents())
        self.assertEqual({"small_deletion"}, self._filtered(node))

    def test_ands_with_the_class_restriction(self):
        """ Orthogonal filters - both restrict, neither joins the damage OR pool """
        node = DamageNode(analysis=self.analysis, variant_class=[VariantClass.DELETION],
                          structural=StructuralFilter.ONLY)
        self.assertEqual({"sv_deletion"}, self._filtered(node))
