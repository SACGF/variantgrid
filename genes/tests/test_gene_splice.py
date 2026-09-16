"""Splice event identity - the alt encoding and how a junction gets its name."""
from django.test import TestCase

from genes.gene_splice import (
    SPLICE_STRING_PATTERN,
    SpliceEventResolver,
    coordinate_label,
    find_splice_events_for_string,
    get_splice_event_variant,
    resolve_splice_string,
)
from genes.models import HGNC, GeneLevelId, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_splice_event_variant
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from snpdb.models import GenomeBuild, Variant

AR_HGNC_ID = 644


class SpliceAltTest(TestCase):
    """ The label segment only a splice alt carries """

    def test_round_trip_with_a_label(self):
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.SPLICE, GeneIdNamespace.HGNC,
                                          AR_HGNC_ID, "V7")
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V7>", alt)
        self.assertEqual((GeneLevelSymbolicAlt.SPLICE, GeneIdNamespace.HGNC, AR_HGNC_ID, "V7"),
                         GeneLevelSymbolicAlt.parse(alt))

    def test_a_coordinate_label_survives_the_round_trip(self):
        """ The label for a junction we have no name for is its coordinates, underscore separated """
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.SPLICE, GeneIdNamespace.HGNC,
                                          AR_HGNC_ID, "X_66905968_66914514")
        self.assertEqual("X_66905968_66914514", GeneLevelSymbolicAlt.parse(alt)[3])

    def test_other_kinds_have_no_label(self):
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.GAIN, GeneIdNamespace.HGNC, AR_HGNC_ID)
        self.assertEqual((GeneLevelSymbolicAlt.GAIN, GeneIdNamespace.HGNC, AR_HGNC_ID, None),
                         GeneLevelSymbolicAlt.parse(alt))

    def test_a_splice_alt_without_a_label_is_not_one(self):
        self.assertIsNone(GeneLevelSymbolicAlt.parse(f"<SPLICE:HGNC:{AR_HGNC_ID}>"))


class SpliceEventResolutionTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=HGNCImport.objects.create(),
                            status=HGNCStatus.APPROVED, approved_name="androgen receptor")

    def test_seeded_junction_resolves_to_its_label(self):
        resolved = SpliceEventResolver(self.genome_build).resolve("AR", "chrX", 66905968, 66914514)
        self.assertEqual("V7", resolved.label)
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V7>", resolved.alt)
        self.assertEqual("AR-V7 splice variant", resolved.display)

    def test_unnamed_junction_is_labelled_with_its_coordinates(self):
        resolved = SpliceEventResolver(self.genome_build).resolve("AR", "chrX", 66905968, 66999999)
        contig = self.genome_build.chrom_contig_mappings["chrX"]
        self.assertEqual(coordinate_label(contig, 66905968, 66999999), resolved.label)
        self.assertEqual("AR X_66905968_66999999", resolved.display,
                         "reads as raw coordinates, which is the prompt to name it")

    def test_two_events_in_one_gene_are_two_variants(self):
        first = create_splice_event_variant("AR", "V7")
        second = create_splice_event_variant("AR", "V9")
        self.assertNotEqual(first.variant, second.variant)
        self.assertEqual(first.variant.locus, second.variant.locus, "same gene, so same locus")

    def test_variant_reads_back_as_its_event(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        read_back = get_splice_event_variant(splice_event_variant.variant)
        self.assertEqual("V7", read_back.label)
        self.assertEqual("AR", read_back.gene.symbol_str)
        self.assertEqual("AR V7", read_back.canonical_str)
        self.assertEqual("AR-V7 splice variant", read_back.display,
                         "the seeded name, found on gene symbol and label")

    def test_a_junction_we_have_no_name_for_displays_its_label(self):
        SpliceEvent.objects.filter(label="V7").delete()
        self.assertEqual("AR V7", create_splice_event_variant("AR", "V7").display)


class SpliceStringResolutionTest(TestCase):
    """ What a lab writes as a classification target, and what search runs on """

    @classmethod
    def setUpTestData(cls):
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(AR_HGNC_ID, "AR"), (3236, "EGFR"), (7029, "MET")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

    def test_every_written_form_reaches_one_alt(self):
        """ Including with the spaces removed, which is how an imported c.HGVS arrives """
        for written in ["AR V7", "ARV7", "ar v7", "AR-V7 splice variant", "AR-V7splicevariant"]:
            self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V7>",
                             resolve_splice_string(written).alt, written)

    def test_the_display_a_report_writes_resolves(self):
        resolved = resolve_splice_string("METexon14skipping")
        self.assertEqual("MET ex14skip", resolved.canonical_str)
        self.assertEqual("MET exon 14 skipping", resolved.display)

    def test_a_label_we_have_no_row_for_mints_nothing(self):
        self.assertIsNone(resolve_splice_string("AR V9"))
        self.assertIsNone(resolve_splice_string("AR-V9 splice variant"))

    def test_a_gene_we_do_not_know_mints_nothing(self):
        self.assertIsNone(resolve_splice_string("NOTAGENE V7"))
        self.assertIsNone(resolve_splice_string("not a splice event"))

    def test_the_coordinate_form_needs_the_variant_to_already_exist(self):
        """ A junction nothing names is something we loaded, not something free text may mint """
        self.assertIsNone(resolve_splice_string("AR X_66905968_66999999"))
        create_splice_event_variant("AR", "X_66905968_66999999")
        resolved = resolve_splice_string("AR X_66905968_66999999")
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:X_66905968_66999999>", resolved.alt)
        self.assertIsNone(resolve_splice_string("ARX_66905968_66999998"),
                          "a junction we did not load")

    def test_resolves_to_the_coordinate_the_loader_writes(self):
        """ A classification target and the loader put the same coordinate through the pipeline, so
            both arrive at one Variant """
        from_loader = SpliceEventResolver(GenomeBuild.grch37()).resolve("AR", "chrX", 66905968, 66914514)
        self.assertEqual(from_loader.variant_coordinate, resolve_splice_string("AR V7").variant_coordinate)

    def test_the_pattern_gates_what_search_runs_on(self):
        for written in ["AR V7", "EGFRvIII", "MET exon 14 skipping", "METex14skip",
                        "AR X_66905968_66914514"]:
            self.assertTrue(SPLICE_STRING_PATTERN.match(written), written)
        for written in ["BCR::ABL1", "EGFR amplification", "NM_000059.3:c.1234A>G"]:
            self.assertFalse(SPLICE_STRING_PATTERN.match(written), written)


class SpliceEventLookupTest(TestCase):
    """ find_splice_events_for_string is what search runs on, so it must never create anything """

    @classmethod
    def setUpTestData(cls):
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=HGNCImport.objects.create(),
                            status=HGNCStatus.APPROVED, approved_name="androgen receptor")

    def test_finds_an_existing_event_by_any_written_form(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        for written in ["AR V7", "ARV7", "AR-V7 splice variant"]:
            found = find_splice_events_for_string(written)
            self.assertEqual([splice_event_variant.variant], [s.variant for s in found], written)

    def test_creates_nothing(self):
        before = (Variant.objects.count(), GeneLevelId.objects.count())
        self.assertEqual([], find_splice_events_for_string("AR V7"))
        self.assertEqual([], find_splice_events_for_string("NOTAGENE V7"))
        self.assertEqual(before, (Variant.objects.count(), GeneLevelId.objects.count()))
