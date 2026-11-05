""" The Effect node's variant type restriction - the All Variants grid's types (#2674) """
from importlib import import_module

from django.test import SimpleTestCase, TestCase

from analysis.models.nodes.filters.damage_node import DamageNode
from analysis.tests.utils import AnalysisSetupMixin, annotate_variants, node_keeps
from annotation.enums import PathogenicityImpact
from library.genomics.vcf_enums import VariantClass, VCFSymbolicAllele
from snpdb.models.models_variant import VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)
from snpdb.variant_filters import VariantType, get_all_variant_types

variant_class_to_variant_types = import_module(
    "analysis.migrations.0014_damagenode_variant_types").variant_class_to_variant_types


class DamageNodeVariantTypesTest(AnalysisSetupMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        sv_deletion = VariantCoordinate(chrom="1", position=10_000, ref="A", alt=VCFSymbolicAllele.DEL,
                                        end=1_010_000, svlen=-1_000_000)
        cls.variants = annotate_variants(cls.analysis.annotation_version.variant_annotation_version, {
            "snv_high": (slowly_create_test_variant("1", 1000, "A", "T", cls.grch37),
                         {"impact": PathogenicityImpact.HIGH}),
            "snv_low": (slowly_create_test_variant("1", 1100, "A", "T", cls.grch37),
                        {"impact": PathogenicityImpact.LOW}),
            "insertion": (slowly_create_test_variant("1", 1200, "A", "AT", cls.grch37), {}),
            "small_deletion": (slowly_create_test_variant("1", 3_000_000, "AT", "A", cls.grch37), {}),
            "sv_deletion": (slowly_create_test_variant_from_coordinate(sv_deletion, cls.grch37), {}),
        })

    def _filtered(self, node: DamageNode) -> set[str]:
        return node_keeps(node, self.variants)

    def test_empty_means_every_type(self):
        node = DamageNode(analysis=self.analysis, variant_types=[])
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_every_type_means_every_type(self):
        node = DamageNode(analysis=self.analysis, variant_types=get_all_variant_types())
        self.assertFalse(node.modifies_parents())
        self.assertIsNone(node._get_node_q())

    def test_restricts_by_stored_alt(self):
        node = DamageNode(analysis=self.analysis, variant_types=[VariantType.SNV])
        self.assertTrue(node.modifies_parents())
        self.assertEqual({"snv_high", "snv_low"}, self._filtered(node))

    def test_small_and_structural_deletions_are_different_types(self):
        node = DamageNode(analysis=self.analysis, variant_types=[VariantType.INDEL])
        self.assertEqual({"insertion", "small_deletion"}, self._filtered(node))
        node = DamageNode(analysis=self.analysis, variant_types=[VCFSymbolicAllele.DEL])
        self.assertEqual({"sv_deletion"}, self._filtered(node))

    def test_ands_with_scoring_pool(self):
        """ It restricts rather than joining the damage OR pool """
        node = DamageNode(analysis=self.analysis, variant_types=[VariantType.SNV], impact_min=PathogenicityImpact.HIGH)
        self.assertEqual({"snv_high"}, self._filtered(node))


class DamageNodeTypeChipsTest(SimpleTestCase):
    @staticmethod
    def _chip_texts(**kwargs) -> list[str]:
        return [chip.text for chip in DamageNode(**kwargs).get_node_chips()]

    def test_no_restriction_no_chips(self):
        self.assertEqual([], self._chip_texts(variant_types=[]))
        self.assertEqual([], self._chip_texts(variant_types=get_all_variant_types()))

    def test_lists_the_types_in_grid_order(self):
        self.assertEqual(["SNV, Indel, DEL"],
                         self._chip_texts(variant_types=[VCFSymbolicAllele.DEL, VariantType.INDEL, VariantType.SNV]))


class VariantClassMigrationTest(SimpleTestCase):
    """ The 0014 data migration's mapping from VEP classes plus the structural radio """
    ALL = [VariantType.SNV, VariantType.INDEL, VariantType.COMPLEX, "<DEL>", "<DUP>",
           VariantType.FUSION, VariantType.COPY_NUMBER, VariantType.SPLICE]

    def _mapped(self, variant_class=(), exclude=False, structural="A") -> list[str]:
        return variant_class_to_variant_types(list(variant_class), exclude, structural, self.ALL)

    def test_nothing_set(self):
        self.assertEqual([], self._mapped())

    def test_classes(self):
        self.assertEqual([VariantType.SNV], self._mapped([VariantClass.SNV]))
        # VEP's deletion covers both the small and the symbolic form
        self.assertEqual([VariantType.INDEL, VariantType.COMPLEX, "<DEL>"], self._mapped([VariantClass.DELETION]))
        # An alt the deployment doesn't offer is dropped
        self.assertEqual([VariantType.INDEL, VariantType.COMPLEX], self._mapped([VariantClass.INSERTION]))
        self.assertEqual(["<DUP>", VariantType.COPY_NUMBER], self._mapped([VariantClass.COPY_NUMBER_GAIN]))
        self.assertEqual(["<DEL>", "<DUP>"], self._mapped([VariantClass.INTERCHROMOSOMAL_TRANSLOCATION]))
        self.assertEqual([VariantType.FUSION], self._mapped([VariantClass.GENE_FUSION]))

    def test_exclude_is_the_complement(self):
        self.assertNotIn(VariantType.SNV, self._mapped([VariantClass.SNV], exclude=True))
        self.assertIn(VariantType.INDEL, self._mapped([VariantClass.SNV], exclude=True))

    def test_structural_radio(self):
        self.assertEqual(["<DEL>", "<DUP>", VariantType.FUSION, VariantType.COPY_NUMBER, VariantType.SPLICE],
                         self._mapped(structural="O"))
        self.assertEqual([VariantType.SNV, VariantType.INDEL, VariantType.COMPLEX], self._mapped(structural="E"))
        # The radio narrows the classes chosen: deletions with structural only were the SV deletions
        self.assertEqual(["<DEL>"], self._mapped([VariantClass.DELETION], structural="O"))
        self.assertEqual([VariantType.INDEL, VariantType.COMPLEX], self._mapped([VariantClass.DELETION], structural="E"))
