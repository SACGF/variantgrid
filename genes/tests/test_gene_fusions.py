"""Gene fusion identity - resolution, the alt encoding, and what collapses onto one Variant."""
from typing import Optional

from django.core.exceptions import ValidationError
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import _create_fake_gene_version, _insert_transcript_data
from genes.gene_fusions import (
    GeneFusionResolver,
    find_gene_fusions_for_string,
    resolve_fusion_string,
)
from genes.models import (
    HGNC,
    GeneLevelId,
    GeneFusion,
    GeneSymbol,
    GeneSymbolAlias,
    HGNCImport,
    ReleaseTranscriptVersion,
)
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource, HGNCStatus
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from snpdb.models import Contig, GenomeBuild, SequenceRole


def _make_release_gene(genome_build, release, gene_id, gene_symbol, transcript_id, contig, start,
                       hgnc_id=None):
    """ A gene of the release with one transcript, so a breakpoint inside it resolves """
    gene_version = _create_fake_gene_version(genome_build, gene_id, gene_symbol,
                                             AnnotationConsortium.ENSEMBL)
    gene_version.hgnc_id = hgnc_id
    gene_version.save()
    data = {
        "id": transcript_id,
        "gene_name": gene_symbol,
        "biotype": [],
        "genome_builds": {
            genome_build.name: {
                "url": "fake",
                "exons": [[start, start + 10_000, 0, 1, 10_001, None]],
                "contig": contig,
                "strand": "+",
                "cds_end": start + 10_000,
                "cds_start": start,
            }
        },
    }
    transcript_version = _insert_transcript_data(genome_build, data, gene_version, release)
    ReleaseTranscriptVersion.objects.get_or_create(release=release, transcript_version=transcript_version)
    return gene_version.gene


class GeneFusionTestCase(TestCase):
    """ The test data's genes, plus the two cases that make identity interesting: SEPT14, which HGNC
        renamed, and the clone-based identifiers that have no HGNC entry at all """

    @classmethod
    def setUpTestData(cls):
        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(1097, "BRAF"), (1014, "BCR"), (76, "ABL1"), (3236, "EGFR"),
                           (10261, "ROS1"), (17127, "GOPC"), (1697, "CD74"), (7989, "NOTCH2NL"),
                           (9236, "PPARG"), (8622, "PAX8"), (3477, "ETV6"), (8033, "NTRK3"),
                           (3446, "ENTPD3"), (10305, "RPL14"), (7527, "MYC"), (9410, "PVT1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk

        # SEPT14 was renamed SEPTIN14 - the file still says the old name
        GeneSymbol.objects.get_or_create(symbol="SEPTIN14")
        HGNC.objects.create(pk=17167, gene_symbol_id="SEPTIN14", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="septin 14",
                            previous_symbols="SEPT14")
        GeneSymbolAlias.objects.create(alias="SEPT14", gene_symbol_id="SEPTIN14",
                                       source=GeneSymbolAliasSource.HGNC)
        cls.hgnc_ids["SEPTIN14"] = 17167

        # ACPP is HGNC's previous symbol for ACP3, and Ensembl still writes it, so there is a
        # GeneSymbol row for it and no GeneSymbolAlias - the matcher stops at the direct hit
        GeneSymbol.objects.get_or_create(symbol="ACP3")
        GeneSymbol.objects.get_or_create(symbol="ACPP")
        HGNC.objects.create(pk=125, gene_symbol_id="ACP3", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="acid phosphatase 3",
                            previous_symbols="ACPP")
        cls.hgnc_ids["ACP3"] = 125

        # SEPT2 is SEPTIN2's previous symbol and, on an unrelated gene, one of SEPTIN6's aliases -
        # GeneSymbolAlias holds a row for each, and keeps only one of them per alias
        for pk, symbol, previous, alias in [(7729, "SEPTIN2", "SEPT2", ""),
                                            (15848, "SEPTIN6", "SEPT6", "SEP2,SEPT2")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name",
                                previous_symbols=previous, alias_symbols=alias)
            cls.hgnc_ids[symbol] = pk
        GeneSymbolAlias.objects.create(alias="SEPT2", gene_symbol_id="SEPTIN6",
                                       source=GeneSymbolAliasSource.HGNC)

        # Clone-based identifiers: known to the annotation as gene symbols, but HGNC carries neither
        for symbol in ["RP11-458D21.5", "AC016683.6"]:
            GeneSymbol.objects.get_or_create(symbol=symbol)

    def setUp(self):
        self.resolver = GeneFusionResolver()


class TestGeneLevelContig(GeneFusionTestCase):

    def test_contig_exists_for_every_build(self):
        contig = Contig.get_gene_level()
        self.assertEqual(SequenceRole.VG_GENE_LEVEL_FAKE_CONTIG, contig.role)
        self.assertTrue(contig.is_gene_level)
        build_names = set(contig.get_genome_builds(require_annotation=False).values_list("name", flat=True))
        self.assertIn("GRCh37", build_names)
        self.assertIn("GRCh38", build_names)

    def test_not_a_standard_contig(self):
        """ standard_contigs is what VCF chrom mappings and header writing use """
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.assertNotIn(Contig.get_gene_level(), list(genome_build.standard_contigs))


class TestFusionResolution(GeneFusionTestCase):

    def test_hgnc_id_is_the_identity(self):
        gene = self.resolver.resolve_side("BRAF")
        self.assertEqual(self.hgnc_ids["BRAF"], gene.gene_level_id.pk)
        self.assertFalse(gene.gene_level_id.is_custom)
        self.assertEqual(GeneIdNamespace.HGNC, gene.gene_level_id.alt_namespace)

    def test_renamed_symbol_resolves_to_approved(self):
        """ The file says SEPT14; HGNC says SEPTIN14 """
        gene = self.resolver.resolve_side("SEPT14")
        self.assertEqual("SEPTIN14", gene.resolved_symbol)
        self.assertEqual(self.hgnc_ids["SEPTIN14"], gene.gene_level_id.pk)
        self.assertTrue(gene.was_renamed)

    def test_multi_gene_side_prefers_the_hgnc_member(self):
        gene = self.resolver.resolve_side("ROS1;GOPC")
        self.assertEqual("ROS1", gene.resolved_symbol)
        self.assertEqual("ROS1;GOPC", gene.written)

    def test_multi_gene_side_skips_the_clone_identifier(self):
        gene = self.resolver.resolve_side("RP11-458D21.5;NOTCH2NL")
        self.assertEqual("NOTCH2NL", gene.resolved_symbol)

    def test_slash_separates_a_multi_gene_side(self):
        gene = self.resolver.resolve_side("PPARG/AC016683.6")
        self.assertEqual("PPARG", gene.resolved_symbol)

    def test_hyphen_does_not_split_a_clone_identifier(self):
        """ RP11-458D21.5 contains a hyphen, so a hyphen can't separate names within a cell """
        self.assertEqual(["RP11-458D21.5"], GeneFusionResolver.split_gene_names("RP11-458D21.5"))

    def test_gene_with_no_hgnc_gets_a_local_id(self):
        gene = self.resolver.resolve_side("RP11-458D21.5")
        gene_level_id = gene.gene_level_id
        self.assertTrue(gene_level_id.is_custom)
        self.assertGreaterEqual(gene_level_id.pk, GeneLevelId.CUSTOM_ID_START)
        self.assertEqual("RP11-458D21.5", gene_level_id.gene_symbol_id)
        self.assertIsNone(gene_level_id.hgnc_id)
        self.assertEqual(GeneIdNamespace.GENE, gene_level_id.alt_namespace)

    def test_unknown_name_still_gets_an_identity(self):
        """ A caller naming only symbols we've never seen still described a real event """
        gene = self.resolver.resolve_side("CTD-2035E11.3")
        self.assertTrue(gene.gene_level_id.is_custom)
        self.assertEqual("CTD-2035E11.3", gene.gene_level_id.symbol_str)
        self.assertIsNone(gene.gene_level_id.gene_symbol_id)

    def test_previous_symbol_resolves_when_a_gene_symbol_row_shadows_it(self):
        """ ACPP is a GeneSymbol in its own right, so the matcher never reaches an alias - HGNC's
            own previous symbols are what places it on ACP3 """
        gene = self.resolver.resolve_side("ACPP")
        self.assertEqual("ACP3", gene.resolved_symbol)
        self.assertEqual(self.hgnc_ids["ACP3"], gene.gene_level_id.pk)
        self.assertTrue(gene.was_renamed)

    def test_a_rename_outranks_another_gene_alias_of_the_same_name(self):
        """ SEPT2 became SEPTIN2; SEPTIN6 merely lists it as an alias, and GeneSymbolAlias has kept
            that one. A rename says 'the same gene', an alias does not """
        gene = self.resolver.resolve_side("SEPT2")
        self.assertEqual("SEPTIN2", gene.resolved_symbol)

    def test_local_ids_are_allocated_without_collision(self):
        first = self.resolver.resolve_side("RP11-458D21.5").gene_level_id
        second = self.resolver.resolve_side("AC016683.6").gene_level_id
        self.assertNotEqual(first.pk, second.pk)
        self.assertEqual(first.pk, self.resolver.resolve_side("RP11-458D21.5").gene_level_id.pk)


class TestFusionVariant(GeneFusionTestCase):

    def _fusion(self, gene_a: Optional[str], gene_b: Optional[str], directionality_known=True) -> GeneFusion:
        return create_gene_fusion(gene_a, gene_b, directionality_known, resolver=self.resolver)

    def test_variant_is_on_the_gene_level_contig(self):
        gene_fusion = self._fusion("BCR", "ABL1")
        variant = gene_fusion.variant
        self.assertTrue(variant.is_gene_level)
        self.assertEqual(Contig.get_gene_level(), variant.locus.contig)
        self.assertEqual(self.hgnc_ids["BCR"], variant.locus.position)
        self.assertEqual("N", variant.locus.ref.seq)
        self.assertEqual(0, variant.svlen, "svlen is 0 not null, so unique_together actually constrains")

    def test_alt_carries_the_partner(self):
        gene_fusion = self._fusion("BCR", "ABL1")
        kind, namespace, partner_id = GeneLevelSymbolicAlt.parse(gene_fusion.variant.alt.seq)
        self.assertEqual(GeneLevelSymbolicAlt.FUSION, kind)
        self.assertEqual(GeneIdNamespace.HGNC, namespace)
        self.assertEqual(self.hgnc_ids["ABL1"], partner_id)

    def test_reciprocal_fusions_are_different_variants(self):
        """ The 5' promoter drives a different protein - BCR-ABL1 is not ABL1-BCR """
        bcr_abl1 = self._fusion("BCR", "ABL1")
        abl1_bcr = self._fusion("ABL1", "BCR")
        self.assertNotEqual(bcr_abl1.variant, abl1_bcr.variant)
        self.assertNotEqual(bcr_abl1, abl1_bcr)

    def test_same_pair_collapses_to_one_fusion(self):
        first = self._fusion("ENTPD3", "RPL14")
        second = self._fusion("ENTPD3", "RPL14")
        self.assertEqual(first, second)
        self.assertEqual(1, GeneFusion.objects.filter(anchor_id=self.hgnc_ids["ENTPD3"]).count())

    def test_unordered_anchors_on_the_lower_id_either_way(self):
        forward = self._fusion("ENTPD3", "RPL14", directionality_known=False)
        reverse = self._fusion("RPL14", "ENTPD3", directionality_known=False)
        self.assertEqual(forward, reverse)
        self.assertFalse(forward.is_ordered)
        self.assertEqual(min(self.hgnc_ids["ENTPD3"], self.hgnc_ids["RPL14"]), forward.anchor_id)

    def test_ordered_and_unordered_are_distinct(self):
        ordered = self._fusion("ENTPD3", "RPL14", directionality_known=True)
        unordered = self._fusion("ENTPD3", "RPL14", directionality_known=False)
        self.assertNotEqual(ordered.variant, unordered.variant)

    def test_unknown_partner(self):
        gene_fusion = self._fusion("EGFR", None)
        self.assertIsNone(gene_fusion.partner_id)
        self.assertEqual("<FUSION:UNKNOWN>", gene_fusion.variant.alt.seq)
        self.assertEqual("EGFR::?", gene_fusion.canonical_str)

    def test_only_the_three_prime_side_known_does_not_assert_direction(self):
        """ Anchoring a 3' partner as if it were 5' would claim something the caller didn't """
        gene_fusion = self._fusion(None, "EGFR")
        self.assertFalse(gene_fusion.is_ordered)
        self.assertEqual("<FUSION_UNORDERED:UNKNOWN>", gene_fusion.variant.alt.seq)

    def test_canonical_str_uses_symbols_not_numbers(self):
        """ What has to leave the system, since local ids mean nothing elsewhere """
        gene_fusion = self._fusion("EGFR", "SEPT14")
        self.assertEqual("EGFR::SEPTIN14", gene_fusion.canonical_str)

    def test_clean_rejects_a_mismatched_alt(self):
        gene_fusion = self._fusion("BCR", "ABL1")
        gene_fusion.is_ordered = False
        with self.assertRaises(ValidationError):
            gene_fusion.clean()


class TestBreakpointResolution(TestCase):
    """ Position rather than spelling - the caller's breakpoint against the release's transcripts """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.release = annotation_version.variant_annotation_version.gene_annotation_release

        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(125, "ACP3"), (3490, "ETV1"), (900, "OVER_A"), (901, "OVER_B")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk

        cls.acp3_gene = _make_release_gene(cls.genome_build, cls.release, "ENSG00000000010", "ACP3",
                                           "ENST00000000010.1", "3", 132_000_000, hgnc_id=125)
        # Two genes over one position: the name written is all there is to choose between them
        cls.over_a_gene = _make_release_gene(cls.genome_build, cls.release, "ENSG00000000011", "OVER_A",
                                             "ENST00000000011.1", "3", 140_000_000, hgnc_id=900)
        cls.over_b_gene = _make_release_gene(cls.genome_build, cls.release, "ENSG00000000012", "OVER_B",
                                             "ENST00000000012.1", "3", 140_000_000, hgnc_id=901)

    def setUp(self):
        self.resolver = GeneFusionResolver()

    def test_breakpoint_beats_the_name_written(self):
        """ A name no release and no HGNC carries, at a position that is unambiguously ACP3 """
        gene = self.resolver.resolve_side("WHATEVER1", breakpoint="chr3:132005000",
                                          genome_build=self.genome_build)
        self.assertEqual("ACP3", gene.resolved_symbol)
        self.assertEqual(self.hgnc_ids["ACP3"], gene.gene_level_id.pk)
        self.assertEqual({self.acp3_gene.pk},
                         set(gene.gene_level_id.genes.values_list("pk", flat=True)))

    def test_chromosome_is_resolved_through_the_build(self):
        """ A caller writes chr3 and a VCF writes 3 - both are the same contig """
        gene = self.resolver.resolve_side("WHATEVER2", breakpoint="3:132005000",
                                          genome_build=self.genome_build)
        self.assertEqual("ACP3", gene.resolved_symbol)

    def test_overlapping_genes_are_decided_by_the_name(self):
        gene = self.resolver.resolve_side("OVER_B", breakpoint="chr3:140005000",
                                          genome_build=self.genome_build)
        self.assertEqual("OVER_B", gene.resolved_symbol)
        self.assertEqual({self.over_b_gene.pk},
                         set(gene.gene_level_id.genes.values_list("pk", flat=True)))

    def test_overlapping_genes_with_no_name_match_fall_back_to_the_name(self):
        """ Picking one of two would be a guess, so the caller's name stands and no gene is recorded """
        gene = self.resolver.resolve_side("WHATEVER3", breakpoint="chr3:140005000",
                                          genome_build=self.genome_build)
        self.assertEqual("WHATEVER3", gene.resolved_symbol)
        self.assertFalse(gene.gene_level_id.genes.exists())

    def test_breakpoint_in_nothing_falls_back_to_the_name(self):
        gene = self.resolver.resolve_side("ACP3", breakpoint="chr3:1000",
                                          genome_build=self.genome_build)
        self.assertEqual(self.hgnc_ids["ACP3"], gene.gene_level_id.pk)
        self.assertFalse(gene.gene_level_id.genes.exists())

    def test_no_build_leaves_resolution_on_the_name(self):
        gene = self.resolver.resolve_side("ACP3", breakpoint="chr3:132005000")
        self.assertEqual(self.hgnc_ids["ACP3"], gene.gene_level_id.pk)
        self.assertFalse(gene.gene_level_id.genes.exists())


class TestFusionString(GeneFusionTestCase):
    """ How a classification target ('BCR::ABL1') reaches the same fusion the loader made """

    def test_separators(self):
        expected = ("BCR", "ABL1")
        for fusion_string in ["BCR::ABL1", "BCR-ABL1", "BCR~ABL1", "BCR/ABL1", "BCR--ABL1"]:
            self.assertEqual(expected, GeneFusionResolver.split_fusion_string(fusion_string), fusion_string)

    def test_resolves_to_the_coordinate_the_loader_writes(self):
        """ A classification target and the loader put the same coordinate through the pipeline, so
            both arrive at one Variant """
        from_loader = create_gene_fusion("BCR", "ABL1", resolver=self.resolver)
        resolved = resolve_fusion_string("BCR::ABL1")
        self.assertEqual(from_loader.variant.coordinate, resolved.variant_coordinate)

    def test_unknown_genes_do_not_mint_a_fusion(self):
        """ Otherwise any hyphenated string would become a gene fusion """
        self.assertIsNone(resolve_fusion_string("SOME-JUNK"))
        self.assertIsNone(resolve_fusion_string("not a fusion at all"))

    def test_resolving_a_string_creates_no_variant(self):
        """ The Variant is the insert pipeline's to make - @see snpdb.gene_level_variants """
        before = GeneFusion.objects.count()
        self.assertIsNotNone(resolve_fusion_string("BCR::ABL1"))
        self.assertEqual(before, GeneFusion.objects.count())


class TestFusionLookup(GeneFusionTestCase):
    """ find_gene_fusions_for_string is what search runs on, so it must never create anything """

    def test_finds_an_existing_fusion(self):
        gene_fusion = create_gene_fusion("BCR", "ABL1", resolver=self.resolver)
        self.assertEqual([gene_fusion], find_gene_fusions_for_string("BCR::ABL1"))

    def test_finds_either_direction(self):
        create_gene_fusion("BCR", "ABL1", resolver=self.resolver)
        self.assertEqual(1, len(find_gene_fusions_for_string("ABL1::BCR")))

    def test_resolves_the_alias(self):
        gene_fusion = create_gene_fusion("EGFR", "SEPTIN14", resolver=self.resolver)
        self.assertEqual([gene_fusion], find_gene_fusions_for_string("EGFR::SEPT14"))

    def test_creates_nothing(self):
        before = (GeneFusion.objects.count(), GeneLevelId.objects.count())
        self.assertEqual([], find_gene_fusions_for_string("BCR::ABL1"))
        self.assertEqual([], find_gene_fusions_for_string("NOTAGENE::ALSONOTAGENE"))
        self.assertEqual(before, (GeneFusion.objects.count(), GeneLevelId.objects.count()))
