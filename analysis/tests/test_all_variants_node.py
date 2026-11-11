""" AllVariantsNode's variant type checkboxes - the Q is built as an OR of the selected types or as the
    complement of the unselected ones, whichever is cheaper, and both have to select the same variants -
    and its chromosome and cytoband pickers """
from importlib import import_module

from django.test import SimpleTestCase, TestCase, override_settings

from analysis.forms.forms_nodes import AllVariantsNodeForm
from analysis.models import AllVariantsNode
from analysis.tests.utils import AnalysisSetupMixin
from annotation.fake_data import create_fake_variant_gene_overlaps
from genes.fake_data import create_fake_transcript_version, create_gene_level_variant
from genes.models import GeneSymbol
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from library.genomics.vcf_enums import GeneLevelSymbolicAlt, VCFSymbolicAllele
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME, GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models import Variant, VariantCoordinate
from snpdb.models.models_genome import Cytoband
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)
from snpdb.variant_filters import GENE_LEVEL_VARIANT_TYPES, VariantType, get_all_variant_types

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
        # A gene-level variant is on a chromosome when one of its genes is - put them all on RUNX1 (chr21)
        cls.gene_contig = cls._add_gene_overlaps(cls.grch37, cls.analysis.annotation_version, [
            v for vt in GENE_LEVEL_VARIANT_TYPES for v in cls.variants_by_type[vt]])

    @staticmethod
    def _add_gene_overlaps(genome_build, annotation_version, variants):
        transcript_version = create_fake_transcript_version(genome_build)
        create_fake_variant_gene_overlaps(annotation_version.variant_annotation_version,
                                          transcript_version.gene_version.gene, variants)
        return transcript_version.contig

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
        """ #2713 - with chromosomes chosen, ticking nothing still keeps the gene-level rows (on the chromosome
            of their gene) and the reference calls, the same as ticking everything """
        chrom_contig = self.variants_by_type[VariantType.SNV][0].locus.contig
        every_variant_id = self._variant_ids(AllVariantsNode.all_variant_types())
        for description, variant_types in {"nothing": [], "everything": AllVariantsNode.all_variant_types()}.items():
            with self.subTest(description):
                node = AllVariantsNode.objects.create(analysis=self.analysis, variant_types=variant_types)
                node.allvariantsnodecontig_set.create(contig=chrom_contig)
                node.allvariantsnodecontig_set.create(contig=self.gene_contig)
                node = AllVariantsNode.objects.get(pk=node.pk)  # contig_ids is cached by the save
                qs = Variant.objects.filter(pk__in=every_variant_id)
                for q in node._get_node_arg_q_dict()[None].values():
                    qs = qs.filter(q)
                self.assertEqual(every_variant_id, set(qs.values_list("pk", flat=True)))


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AllVariantsNodeCytobandTest(AnalysisSetupMixin, TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        contig_11 = cls.grch37.chrom_contig_mappings["11"]
        cls.band_11p13 = Cytoband.objects.get(contig=contig_11, name="p13")  # 31000001-36400000
        cls.in_band = slowly_create_test_variant("11", 31000001, "A", "G", cls.grch37)
        cls.outside_band = slowly_create_test_variant("11", 36400001, "A", "G", cls.grch37)
        cls.on_chrom_2 = slowly_create_test_variant("2", 1000, "A", "G", cls.grch37)
        for symbol in ["CD74", "ROS1"]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
        cls.fusion = create_gene_fusion("CD74", "ROS1").variant
        cls.gene_contig = AllVariantsNodeVariantTypeTest._add_gene_overlaps(
            cls.grch37, cls.analysis.annotation_version, [cls.fusion])
        # The fake RUNX1 transcript is at 21:34787801-35049344
        cls.band_with_gene = Cytoband.objects.get(contig=cls.gene_contig, name="q22.11")

    def _save_form(self, contigs=(), cytobands=()) -> AllVariantsNode:
        node = AllVariantsNode.objects.create(analysis=self.analysis)
        data = {"contigs": [c.pk for c in contigs], "cytobands": [c.pk for c in cytobands],
                "variant_types": get_all_variant_types(), "max_variant": self.on_chrom_2.pk,
                "min_het_or_hom_count": 0}
        form = AllVariantsNodeForm(data, instance=node, genome_build=self.grch37)
        self.assertTrue(form.is_valid(), form.errors)
        form.save()
        return AllVariantsNode.objects.get(pk=node.pk)

    def _filtered_variants(self, node: AllVariantsNode) -> set[Variant]:
        test_variants = [self.in_band, self.outside_band, self.on_chrom_2, self.fusion]
        qs = Variant.objects.filter(pk__in=[v.pk for v in test_variants])
        for q in node._get_node_arg_q_dict()[None].values():
            qs = qs.filter(q)
        return set(qs)

    def test_band_filters_to_its_coordinates(self):
        node = self._save_form(cytobands=[self.band_11p13])
        self.assertEqual(self._filtered_variants(node), {self.in_band})
        self.assertEqual(node.get_node_name(), "11p13")

    def test_chromosomes_or_bands(self):
        contig_2 = self.grch37.chrom_contig_mappings["2"]
        node = self._save_form(contigs=[contig_2], cytobands=[self.band_11p13])
        self.assertEqual(self._filtered_variants(node), {self.in_band, self.on_chrom_2})
        self.assertEqual(node.get_node_name(), "2, 11p13")

    def test_fusion_in_the_band_or_chromosome_of_its_gene(self):
        """ #2769 - a gene-level variant has no coordinate, so it's placed by its genes """
        for description, kwargs in {"band": {"cytobands": [self.band_with_gene]},
                                    "chromosome": {"contigs": [self.gene_contig]}}.items():
            with self.subTest(description):
                self.assertIn(self.fusion, self._filtered_variants(self._save_form(**kwargs)))


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
