""" AllVariantsNode's variant type checkboxes - the Q is built as an OR of the selected types or as the
    complement of the unselected ones, whichever is cheaper, and both have to select the same variants """
from importlib import import_module

from django.test import SimpleTestCase, TestCase, override_settings

from analysis.models import AllVariantsNode
from analysis.tests.utils import AnalysisSetupMixin
from genes.fake_data import create_gene_level_variant
from genes.models import GeneSymbol
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from library.genomics.vcf_enums import GeneLevelSymbolicAlt, VCFSymbolicAllele
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME, GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models import Variant, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)
from snpdb.variant_filters import VariantType, get_all_variant_types

booleans_to_variant_types = import_module(
    "analysis.migrations.0015_allvariantsnode_variant_types").booleans_to_variant_types


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AllVariantsNodeVariantTypeTest(AnalysisSetupMixin, TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        def create(ref, alt, position):
            return slowly_create_test_variant("1", position, ref, alt, cls.grch37)

        def create_gene_level(kind, label=None):
            alt = GeneLevelSymbolicAlt.format(kind, "HGNC", 1, label=label)
            return create_gene_level_variant(VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=1,
                                                               ref=GENE_LEVEL_REF, alt=alt, svlen=GENE_LEVEL_SVLEN))

        for symbol in ["CD74", "ROS1"]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
        deletion = VariantCoordinate(chrom="1", position=5000, ref="A", alt=VCFSymbolicAllele.DEL, svlen=-1000)
        cls.variants_by_type = {
            VariantType.REFERENCE: [create("A", Variant.REFERENCE_ALT, 1000)],
            VariantType.SNV: [create("A", "G", 1000)],
            VariantType.INDEL: [create("A", "AT", 1000), create("AT", "A", 2000)],
            VariantType.COMPLEX: [create("AT", "GC", 3000)],
            VCFSymbolicAllele.DEL: [slowly_create_test_variant_from_coordinate(deletion, cls.grch37)],
            VariantType.FUSION: [create_gene_fusion("CD74", "ROS1").variant],
            VariantType.COPY_NUMBER: [create_gene_level(GeneLevelSymbolicAlt.GAIN)],
            VariantType.SPLICE: [create_gene_level(GeneLevelSymbolicAlt.SPLICE, label="V_7")],
        }

    def _node_variant_ids(self, variant_types) -> set[int]:
        node = AllVariantsNode(analysis=self.analysis, variant_types=list(variant_types))
        qs = Variant.objects.filter(pk__in=self._variant_ids(self.variants_by_type))
        if q := node._get_node_arg_q_dict()[None].get("variant_types"):
            qs = qs.filter(q)
        return set(qs.values_list("pk", flat=True))

    def _variant_ids(self, variant_types) -> set[int]:
        return {v.pk for vt in variant_types for v in self.variants_by_type.get(vt, [])}

    def test_variant_type_selections(self):
        all_types = AllVariantsNode.all_variant_types()
        default_types = get_all_variant_types()
        selections = {"default": default_types}
        for variant_type in all_types:
            selections[f"only {variant_type}"] = [variant_type]
            selections[f"all but {variant_type}"] = [vt for vt in all_types if vt != variant_type]
        for variant_type in default_types:
            selections[f"default but {variant_type}"] = [vt for vt in default_types if vt != variant_type]

        for description, variant_types in selections.items():
            with self.subTest(description):
                self.assertEqual(self._variant_ids(variant_types), self._node_variant_ids(variant_types))

    def test_nothing_ticked_on_a_chromosome_is_no_restriction(self):
        """ #2713 - with a chromosome chosen, ticking nothing still keeps the gene-level rows (they live on
            their own contig) and the reference calls, the same as ticking everything """
        chrom_contig = self.variants_by_type[VariantType.SNV][0].locus.contig
        every_variant_id = self._variant_ids(AllVariantsNode.all_variant_types())
        for description, variant_types in {"nothing": [], "everything": AllVariantsNode.all_variant_types()}.items():
            with self.subTest(description):
                node = AllVariantsNode.objects.create(analysis=self.analysis, variant_types=variant_types)
                node.allvariantsnodecontig_set.create(contig=chrom_contig)
                node = AllVariantsNode.objects.get(pk=node.pk)  # contig_ids is cached by the save
                self.assertEqual([chrom_contig.pk], node.contig_ids)
                qs = Variant.objects.filter(pk__in=every_variant_id)
                for q in node._get_node_arg_q_dict()[None].values():
                    qs = qs.filter(q)
                self.assertEqual(every_variant_id, set(qs.values_list("pk", flat=True)))


class AllVariantsNodeNameTest(SimpleTestCase):
    def test_default_types_say_nothing(self):
        self.assertEqual("", AllVariantsNode().get_node_name())

    def test_restricted_types_are_named_in_grid_order(self):
        node = AllVariantsNode(variant_types=[VCFSymbolicAllele.DEL, VariantType.SNV, VariantType.REFERENCE])
        self.assertEqual("(type: Reference, SNV, DEL)", node.get_node_name())

    def test_nothing_ticked_says_nothing(self):
        self.assertEqual("", AllVariantsNode(variant_types=[]).get_node_name())


class AllVariantsNodeMigrationTest(SimpleTestCase):
    """ The eight checkboxes an AllVariantsNode had become the grid's types (migration 0015) """
    ALL = [VariantType.REFERENCE, VariantType.SNV, VariantType.INDEL, VariantType.COMPLEX,
           VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP, VariantType.FUSION, VariantType.COPY_NUMBER,
           VariantType.SPLICE]
    FIELDS = ["reference", "snps", "indels", "complex_subsitution", "structural_variants",
              "fusions", "copy_number", "splicing"]

    def _mapped(self, *ticked_fields) -> list[str]:
        ticked = {field: field in ticked_fields for field in self.FIELDS}
        with override_settings(VARIANT_SYMBOLIC_ALT_ENABLED=True,
                               VARIANT_SYMBOLIC_ALT_VALID_TYPES={VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP}):
            return booleans_to_variant_types(ticked, self.ALL)

    def test_default_node_is_everything_but_reference(self):
        self.assertEqual([vt for vt in self.ALL if vt != VariantType.REFERENCE], self._mapped(*self.FIELDS[1:]))

    def test_structural_variants_is_every_symbolic_alt(self):
        self.assertEqual([VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP], self._mapped("structural_variants"))

    def test_nothing_ticked_stays_nothing(self):
        self.assertEqual([], self._mapped())
