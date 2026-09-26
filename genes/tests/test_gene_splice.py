"""Splice event identity - the canonical label, the alt encoding, and how a junction gets its name."""
from django.test import TestCase

from annotation.fake_data import get_fake_annotation_version
from genes.gene_splice import (
    SPLICE_STRING_PATTERN,
    SpliceEventResolver,
    canonical_splice_label,
    coordinate_label,
    display_splice_label,
    find_splice_events_for_string,
    get_splice_event_variant,
    resolve_splice_string,
)
from genes.models import HGNC, GeneLevelId, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from genes.tests.gene_level_test_utils import create_splice_event_variant, make_release_gene
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from snpdb.models import GenomeBuild, Variant, VariantCoordinate

AR_HGNC_ID = 644
EGFR_HGNC_ID = 3236
MET_HGNC_ID = 7029
NOTCH2_HGNC_ID = 7882


class SpliceAltTest(TestCase):
    """ The label segment only a splice alt carries """

    def test_round_trip_with_a_label(self):
        """ The alt is a Sequence, so it is stored upper-cased; the label comes back canonical """
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.SPLICE, GeneIdNamespace.HGNC,
                                          AR_HGNC_ID, "v_7")
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", alt)
        self.assertEqual((GeneLevelSymbolicAlt.SPLICE, GeneIdNamespace.HGNC, AR_HGNC_ID, "v_7"),
                         GeneLevelSymbolicAlt.parse(alt))

    def test_other_kinds_have_no_label(self):
        alt = GeneLevelSymbolicAlt.format(GeneLevelSymbolicAlt.GAIN, GeneIdNamespace.HGNC, AR_HGNC_ID)
        self.assertEqual((GeneLevelSymbolicAlt.GAIN, GeneIdNamespace.HGNC, AR_HGNC_ID, None),
                         GeneLevelSymbolicAlt.parse(alt))

    def test_a_splice_alt_without_a_label_is_not_one(self):
        self.assertIsNone(GeneLevelSymbolicAlt.parse(f"<SPLICE:HGNC:{AR_HGNC_ID}>"))

    def test_a_coordinate_string_in_display_case_is_the_stored_coordinate(self):
        """ A gene-level coordinate string is the same coordinate whatever case it arrives in, so a
            record written before the alt was upper-cased re-matches to the Variant we have """
        variant_coordinate = VariantCoordinate.from_string("GENE_LEVEL:3236-3236 <SPLICE:HGNC:3236:v_iii>",
                                                           None)
        self.assertEqual("<SPLICE:HGNC:3236:V_III>", variant_coordinate.alt)


class CanonicalSpliceLabelTest(TestCase):
    """ Every way of writing one junction's name is one label, so it is one Variant """

    WRITTEN_FORMS = {
        "v_7": ["V7", "-V7", "v 7", "V7 splice variant"],
        "v_iii": ["vIII", "-vIII", "VIII", "vIII splice variant"],
        "v_iva": ["vIVa", "VIVa"],
        "exon_14_skipping": ["exon 14 skipping", "ex14skip", "Exon14Skipping", "ex14 skipping"],
    }

    def test_every_written_form_gives_one_label(self):
        for canonical, written_forms in self.WRITTEN_FORMS.items():
            for written in written_forms:
                self.assertEqual(canonical, canonical_splice_label(written), written)

    def test_a_canonical_label_is_left_alone(self):
        for canonical in [*self.WRITTEN_FORMS, "grch37_x_66905968_66914514"]:
            self.assertEqual(canonical, canonical_splice_label(canonical), canonical)

    def test_breakpoints_are_labelled_under_their_build(self):
        """ The same numbers are different junctions in each build, and a gene-level Variant sits on
            the contig they share - so the build is part of the label """
        for written in ["X_66905968_66914514", "chrX:66905968-66914514"]:
            self.assertEqual("grch37_x_66905968_66914514",
                             canonical_splice_label(written, GenomeBuild.grch37()), written)
        self.assertEqual("grch38_x_66905968_66914514",
                         canonical_splice_label("X_66905968_66914514", GenomeBuild.grch38()))
        self.assertIsNone(canonical_splice_label("X_66905968_66914514"),
                          "breakpoints without a build name no one junction")

    def test_a_shape_we_do_not_accept(self):
        for written in ["", "exon 14", "the second one", "NM_000059.3:c.1234A>G"]:
            self.assertIsNone(canonical_splice_label(written), written)


class DisplaySpliceLabelTest(TestCase):
    """ The canonical label written the way the literature writes it """

    @classmethod
    def setUpTestData(cls):
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=HGNCImport.objects.create(),
                            status=HGNCStatus.APPROVED, approved_name="androgen receptor")
        cls.gene = GeneLevelId.objects.create(symbol_str="AR", hgnc_id=AR_HGNC_ID)

    def test_each_shape_displays(self):
        for label, display in [("v_7", "AR-V7 splice"),
                               ("v_iii", "ARvIII splice"),
                               ("v_iva", "ARvIVa splice"),
                               ("exon_14_skipping", "AR exon 14 skipping"),
                               ("grch37_x_66905968_66914514", "AR GRCh37 X:66905968-66914514 splice")]:
            self.assertEqual(display, display_splice_label(self.gene, label))

    def test_a_display_resolves_back_to_its_label(self):
        """ What a report writes, and what we display, is what a lab writes back as a classification target """
        for label, written in [("v_7", "-V7"), ("v_7", "-V7 splice"), ("v_iii", "vIII"),
                               ("v_iva", "vIVa"), ("v_iva", "vIVa splice"),
                               ("exon_14_skipping", " exon 14 skipping")]:
            self.assertEqual(label, canonical_splice_label(written), written)


class SpliceEventResolutionTest(TestCase):
    """ What the TSO 500 importer resolves a caller's breakpoints to """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(AR_HGNC_ID, "AR"), (MET_HGNC_ID, "MET")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

    def test_seeded_junction_resolves_to_its_label(self):
        resolved = SpliceEventResolver(self.genome_build).resolve("AR", "chrX", 66905968, 66914514)
        self.assertEqual("v_7", resolved.label)
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", resolved.alt)
        self.assertEqual("AR-V7 splice variant", resolved.display)

    def test_unnamed_junction_is_labelled_with_its_breakpoints(self):
        resolved = SpliceEventResolver(self.genome_build).resolve("AR", "chrX", 66905968, 66999999)
        contig = self.genome_build.chrom_contig_mappings["chrX"]
        self.assertEqual(coordinate_label(self.genome_build, contig, 66905968, 66999999), resolved.label)
        self.assertEqual("AR GRCh37 X:66905968-66999999 splice", resolved.display,
                         "reads as raw coordinates, which is the prompt to name it")

    def test_two_events_in_one_gene_are_two_variants(self):
        first = create_splice_event_variant("AR", "V7")
        second = create_splice_event_variant("AR", "V9")
        self.assertNotEqual(first.variant, second.variant)
        self.assertEqual(first.variant.locus, second.variant.locus, "same gene, so same locus")

    def test_variant_reads_back_as_its_event(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        read_back = get_splice_event_variant(splice_event_variant.variant)
        self.assertEqual("v_7", read_back.label)
        self.assertEqual("AR", read_back.gene.symbol_str)
        self.assertEqual("AR-V7 splice", read_back.canonical_str)
        self.assertEqual("AR-V7 splice variant", read_back.display,
                         "the seeded wording, found on gene symbol and label")

    def test_a_variant_stored_upper_case_displays_the_canonical_label(self):
        """ The alt holds MET's label as EXON_14_SKIPPING; the report writes the row's name (#1835) """
        splice_event_variant = create_splice_event_variant("MET", "ex14skip")
        self.assertEqual(f"<SPLICE:HGNC:{MET_HGNC_ID}:EXON_14_SKIPPING>",
                         splice_event_variant.variant.alt.seq)
        read_back = get_splice_event_variant(splice_event_variant.variant)
        self.assertEqual("MET exon 14 skipping", read_back.display)
        self.assertEqual("MET exon 14 skipping", read_back.canonical_str)

    def test_a_junction_we_have_no_name_for_displays_its_label(self):
        SpliceEvent.objects.filter(label="v_7").delete()
        self.assertEqual("AR-V7 splice", create_splice_event_variant("AR", "V7").display)


class SpliceJunctionLocusTest(TestCase):
    """ The coordinates an IGV link off the variant page is built from (#1908) """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        cls.other_build = GenomeBuild.grch38()
        hgnc_import = HGNCImport.objects.create()
        GeneSymbol.objects.get_or_create(symbol="AR")
        HGNC.objects.create(pk=AR_HGNC_ID, gene_symbol_id="AR", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="AR approved name")

    def test_a_named_junction_takes_the_builds_splice_event_row(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        self.assertEqual("X:66905968-66914514",
                         splice_event_variant.junction_locus(self.genome_build))

    def test_a_named_junction_has_a_locus_per_build(self):
        splice_event_variant = create_splice_event_variant("AR", "V7")
        self.assertNotEqual(splice_event_variant.junction_locus(self.genome_build),
                            splice_event_variant.junction_locus(self.other_build),
                            "the same junction is different coordinates in each build")

    def test_a_junction_named_by_its_breakpoints_reads_them_off_its_label(self):
        contig = self.genome_build.chrom_contig_mappings["chrX"]
        label = coordinate_label(self.genome_build, contig, 66905968, 66999999)
        splice_event_variant = create_splice_event_variant("AR", label)
        self.assertEqual("X:66905968-66999999",
                         splice_event_variant.junction_locus(self.genome_build))
        self.assertIsNone(splice_event_variant.junction_locus(self.other_build),
                          "breakpoints mean nothing in a build they weren't called in")

    def test_a_junction_with_no_coordinates_anywhere(self):
        SpliceEvent.objects.filter(label="v_7").delete()
        self.assertIsNone(create_splice_event_variant("AR", "V7").junction_locus(self.genome_build))


class SpliceJunctionResolutionTest(TestCase):
    """ A junction whose caller names no gene - SpliceGirl's <DEL> from POS to END """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(AR_HGNC_ID, "AR"), (NOTCH2_HGNC_ID, "NOTCH2"), (900, "OVER_A"), (901, "OVER_B")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

        annotation_version = get_fake_annotation_version(cls.genome_build)
        release = annotation_version.variant_annotation_version.gene_annotation_release
        cls.notch2_gene = make_release_gene(cls.genome_build, release, "ENSG00000134250", "NOTCH2",
                                            "ENST00000256646.1", "1", 120_460_000, hgnc_id=NOTCH2_HGNC_ID)
        # Overlapping genes: 150,005,000-150,010,000 is in both, 150,010,001-150,015,000 only in OVER_B
        make_release_gene(cls.genome_build, release, "ENSG00000000011", "OVER_A",
                          "ENST00000000011.1", "1", 150_000_000, hgnc_id=900)
        make_release_gene(cls.genome_build, release, "ENSG00000000012", "OVER_B",
                          "ENST00000000012.1", "1", 150_005_000, hgnc_id=901)

    def setUp(self):
        self.resolver = SpliceEventResolver(self.genome_build)

    def test_seeded_junction_needs_no_gene_name(self):
        resolved = self.resolver.resolve_junction("chrX", 66905968, 66914514)
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>", resolved.alt)

    def test_breakpoints_are_taken_in_genomic_order(self):
        self.assertEqual(self.resolver.resolve_junction("chrX", 66905968, 66914514),
                         self.resolver.resolve_junction("chrX", 66914514, 66905968))

    def test_unnamed_junction_takes_the_gene_at_its_donor(self):
        resolved = self.resolver.resolve_junction("chr1", 120464432, 120465258)
        contig = self.genome_build.chrom_contig_mappings["chr1"]
        self.assertEqual(NOTCH2_HGNC_ID, resolved.gene.pk)
        self.assertEqual(coordinate_label(self.genome_build, contig, 120464432, 120465258), resolved.label)
        self.assertEqual({self.notch2_gene.pk}, set(resolved.gene.genes.values_list("pk", flat=True)))

    def test_no_gene_at_the_junction_is_unresolved(self):
        self.assertIsNone(self.resolver.resolve_junction("chr1", 1000, 2000))

    def test_overlapping_genes_are_decided_by_the_acceptor(self):
        resolved = self.resolver.resolve_junction("chr1", 150_006_000, 150_012_000)
        self.assertEqual(901, resolved.gene.pk)
        self.assertIsNone(self.resolver.resolve_junction("chr1", 150_006_000, 150_008_000),
                          "both genes span the whole junction - picking one would be a guess")

    def test_a_del_expanded_to_sequence_resolves_like_the_symbolic_one(self):
        """ VG3 wrote the sub-1 kb <DEL>s out as explicit sequence, so the END comes off the ref """
        symbolic = VariantCoordinate(chrom="X", position=66905968, ref="T", alt="<DEL>", svlen=-8546)
        explicit = VariantCoordinate(chrom="X", position=66905968, ref="T" * 8547, alt="T")
        expected = self.resolver.resolve_junction("X", 66905968, 66914514)
        self.assertEqual(expected, self.resolver.resolve_variant_coordinate(symbolic))
        self.assertEqual(expected, self.resolver.resolve_variant_coordinate(explicit))


class SpliceStringResolutionTest(TestCase):
    """ What a lab writes as a classification target, and what search runs on """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(AR_HGNC_ID, "AR"), (EGFR_HGNC_ID, "EGFR"), (MET_HGNC_ID, "MET")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

    def test_every_written_form_reaches_one_alt(self):
        """ Including with the spaces removed, as the classifications imported before #1875 hold it """
        for written in ["AR V7", "ARV7", "ar v7", "AR-V7 splice variant", "AR-V7splicevariant"]:
            self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:V_7>",
                             resolve_splice_string(written).resolved.alt, written)

    def test_a_name_nobody_registered_still_mints_its_variant(self):
        """ The lab's label is the identity - no row has to exist first (#1835) """
        resolved = resolve_splice_string("EGFRvIVa").resolved
        self.assertEqual(f"<SPLICE:HGNC:{EGFR_HGNC_ID}:V_IVA>", resolved.alt)
        self.assertIsNone(resolved.splice_event)
        self.assertEqual("EGFRvIVa splice", resolved.display)

    def test_the_report_and_the_caller_reach_one_alt(self):
        """ A classification target and the loader put the same coordinate through the pipeline """
        from_loader = SpliceEventResolver(self.genome_build).resolve("EGFR", "chr7", 55087058, 55223522)
        self.assertEqual(from_loader.variant_coordinate,
                         resolve_splice_string("EGFRvIII").resolved.variant_coordinate)

    def test_the_display_a_report_writes_resolves(self):
        resolved = resolve_splice_string("METexon14skipping").resolved
        self.assertEqual("exon_14_skipping", resolved.label)
        self.assertEqual("MET exon 14 skipping", resolved.canonical_str)

    def test_a_gene_we_do_not_know_is_refused_with_a_reason(self):
        resolution = resolve_splice_string("NOTAGENE vIII")
        self.assertFalse(resolution)
        self.assertEqual("gene 'NOTAGENE' is not a symbol we know", resolution.reason)

    def test_a_value_that_names_no_junction_is_left_for_the_next_resolver(self):
        for written in ["not a splice event", "NM_000059.3:c.1234A>G"]:
            resolution = resolve_splice_string(written)
            self.assertFalse(resolution.recognised, written)

    def test_breakpoints_resolve_under_the_imported_build(self):
        resolved = resolve_splice_string("AR X_66905968_66999999", GenomeBuild.grch37()).resolved
        self.assertEqual(f"<SPLICE:HGNC:{AR_HGNC_ID}:GRCH37_X_66905968_66999999>", resolved.alt)
        self.assertNotEqual(resolved.alt,
                            resolve_splice_string("AR X_66905968_66999999", GenomeBuild.grch38()).resolved.alt,
                            "the same numbers are a different junction in another build")

    def test_breakpoints_on_a_contig_the_build_does_not_have(self):
        resolution = resolve_splice_string("AR 24_1_2", GenomeBuild.grch37())
        self.assertFalse(resolution)
        self.assertIn("is not a contig", resolution.reason)

    def test_the_pattern_gates_what_search_runs_on(self):
        for written in ["AR V7", "EGFRvIII", "EGFRvIVa", "MET exon 14 skipping", "METex14skip",
                        "AR X_66905968_66914514", "AR chrX:66905968-66914514"]:
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
