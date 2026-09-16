"""Splice event identity - the alt encoding and how a junction gets its name."""
from django.test import TestCase

from genes.gene_splice import SpliceEventResolver, coordinate_label, get_splice_event_variant
from genes.models import HGNC, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_splice_event_variant
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from snpdb.models import GenomeBuild

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
