""" AllVariantsNode's variant type checkboxes - the Q is built as an OR of the selected types or as the
    complement of the unselected ones, whichever is cheaper, and both have to select the same variants """
from django.test import TestCase, override_settings

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
from snpdb.variant_filters import VariantType

NODE_FIELD_FOR_TYPE = {
    VariantType.REFERENCE: "reference",
    VariantType.SNV: "snps",
    VariantType.INDEL: "indels",
    VariantType.COMPLEX: "complex_subsitution",
    VariantType.SYMBOLIC: "structural_variants",
    VariantType.FUSION: "fusions",
    VariantType.COPY_NUMBER: "copy_number",
    VariantType.SPLICE: "splicing",
}


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
            VariantType.SYMBOLIC: [slowly_create_test_variant_from_coordinate(deletion, cls.grch37)],
            VariantType.FUSION: [create_gene_fusion("CD74", "ROS1").variant],
            VariantType.COPY_NUMBER: [create_gene_level(GeneLevelSymbolicAlt.GAIN)],
            VariantType.SPLICE: [create_gene_level(GeneLevelSymbolicAlt.SPLICE, label="V_7")],
        }

    def _node_variant_ids(self, variant_types) -> set[int]:
        node = AllVariantsNode(analysis=self.analysis,
                               **{field: vt in variant_types for vt, field in NODE_FIELD_FOR_TYPE.items()})
        qs = Variant.objects.filter(pk__in=self._variant_ids(AllVariantsNode.VARIANT_TYPES))
        if q := node._get_node_arg_q_dict()[None].get("variant_types"):
            qs = qs.filter(q)
        return set(qs.values_list("pk", flat=True))

    def _variant_ids(self, variant_types) -> set[int]:
        return {v.pk for vt in variant_types for v in self.variants_by_type[vt]}

    def test_variant_type_selections(self):
        all_types = AllVariantsNode.VARIANT_TYPES
        default_types = [vt for vt in all_types if vt != VariantType.REFERENCE]
        selections = {"default": default_types}
        for variant_type in all_types:
            selections[f"only {variant_type}"] = [variant_type]
            selections[f"all but {variant_type}"] = [vt for vt in all_types if vt != variant_type]
        for variant_type in default_types:
            selections[f"default but {variant_type}"] = [vt for vt in default_types if vt != variant_type]

        for description, variant_types in selections.items():
            with self.subTest(description):
                self.assertEqual(self._variant_ids(variant_types), self._node_variant_ids(variant_types))
