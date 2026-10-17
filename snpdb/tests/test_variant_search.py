from unittest import mock

from django.test import SimpleTestCase, TestCase, override_settings

from annotation.cosmic import CosmicAPI
from genes.hgvs import HGVSMatcher
from library.enums.log_level import LogLevel
from library.preview_request import PreviewData
from snpdb.models import GenomeBuild, Sample, Sequence, Variant, VariantCoordinate
from snpdb.search import (
    SearchInput,
    SearchMessage,
    SearchResponse,
    SearchResponsesCombined,
    SearchResult,
    SearchResultMatchStrength,
)
from snpdb.signals.variant_search import (
    _alt_description,
    _alt_mismatch_message,
    _contig_not_in_build_message,
    _gene_symbol_search_message,
    _hgvs_nomenclature_error_message,
    _search_hgvs_using_gene_symbol,
    _transcript_substitution_messages,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class VariantSearchAltDescriptionTestCase(TestCase):
    """ Regression tests for _alt_description (issue #3864).

        When a variant search has no direct hit it looks for "alt alts" - existing Variants at the same
        locus with a different alt - and describes each one. A Variant's alt is a Sequence (not a str), so
        describing it via Sequence.abbreviate(variant.alt) sliced the Sequence object and raised
        "'Sequence' object is not subscriptable" once the alt exceeded the abbreviation length. """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.long_alt = "G" + "A" * 30  # >20 chars, so abbreviate() slices - this is what crashed
        cls.long_alt_variant = slowly_create_test_variant("3", 128198995, "G", cls.long_alt, cls.grch37)
        cls.short_alt_variant = slowly_create_test_variant("3", 128198980, "A", "T", cls.grch37)

    def test_long_alt_on_variant_does_not_crash(self):
        """ Variant.alt is a Sequence - long alts must be abbreviated without subscripting the Sequence. """
        # alt is a Sequence, not a str - this is the condition that triggered the bug
        self.assertIsInstance(self.long_alt_variant.alt, Sequence)
        desc = _alt_description(self.long_alt_variant)
        self.assertIsInstance(desc, str)
        self.assertEqual(desc, Sequence.abbreviate(self.long_alt))

    def test_short_alt_on_variant(self):
        desc = _alt_description(self.short_alt_variant)
        self.assertEqual(desc, "T")

    def test_alt_on_variant_coordinate(self):
        """ VariantCoordinate.alt is a plain str - long alts are still abbreviated. """
        vc = VariantCoordinate(chrom="3", position=128198995, ref="G", alt=self.long_alt)
        self.assertEqual(_alt_description(vc), Sequence.abbreviate(self.long_alt))


class VariantSearchAltMismatchMessageTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.found_variant = slowly_create_test_variant("3", 128198980, "A", "T", cls.grch37)
        cls.searched = VariantCoordinate(chrom="3", position=128198980, ref="A", alt="G")

    def test_hgvs_search_names_both_as_g_hgvs(self):
        message = _alt_mismatch_message(self.searched, self.found_variant, HGVSMatcher.instance(self.grch37))
        self.assertEqual(message, 'No results for "NC_000003.11:g.128198980A>G", '
                                  'but found "NC_000003.11:g.128198980A>T" at the same position')

    def test_coordinate_search_describes_alts(self):
        message = _alt_mismatch_message(self.searched, self.found_variant, None)
        self.assertEqual(message, 'No results for alt "G", but found this using alt "T"')


class CosmicSearchGenomeBuildTestCase(TestCase):
    """ The clinicaltables COSMIC API only has GRCh37/38 data - other builds 500 out (issue #1783) """

    def test_t2t_is_not_supported(self):
        t2t = GenomeBuild.t2tv2()
        self.assertFalse(CosmicAPI.supports_genome_build(t2t))
        with self.assertRaises(ValueError):
            CosmicAPI("COSV53567516", t2t)

    def test_grch37_and_38_are_supported(self):
        for build_name in ["GRCh37", "GRCh38"]:
            self.assertTrue(CosmicAPI.supports_genome_build(GenomeBuild.get_name_or_alias(build_name)))


class SearchAutoOpenGenomeBuildTestCase(TestCase):
    """ A single result outside the preferred build auto-opens only if it can never be in that build (#2174) """

    @classmethod
    def setUpTestData(cls):
        cls.search_input = SearchInput(user=None, search_string="x",
                                       genome_build_preferred=GenomeBuild.get_name_or_alias("GRCh38"))
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _single_preferred_result(self, search_type):
        preview = PreviewData(category=search_type.preview_category(), identifier="x", internal_url="/result",
                              genome_builds={self.grch37}, can_liftover=search_type.preview_can_liftover())
        result = SearchResult(preview=preview, match_strength=SearchResultMatchStrength.STRONG_MATCH)
        response = SearchResponse(search_input=self.search_input, search_type=search_type, matched_pattern=True,
                                  results=[result], total_count=1)
        return SearchResponsesCombined(self.search_input, [response]).single_preferred_result()

    def test_fixed_build_result_in_other_build_auto_opens(self):
        self.assertIsNotNone(self._single_preferred_result(Sample))

    def test_variant_in_other_build_does_not_auto_open(self):
        self.assertIsNone(self._single_preferred_result(Variant))


class HGVSSearchMessagesTestCase(SimpleTestCase):

    def test_transcript_version_used_as_given(self):
        self.assertEqual(_transcript_substitution_messages("NM_000546.6(TP53):c.215C>G", "NM_000546.6"), [])

    def test_transcript_version_added(self):
        messages = _transcript_substitution_messages("NM_000546:c.215C>G", "NM_000546.6")
        self.assertEqual(messages, [SearchMessage('Used transcript version "NM_000546.6"', severity=LogLevel.INFO)])

    def test_transcript_version_substituted(self):
        messages = _transcript_substitution_messages("NM_000546.12(TP53):c.215C>G", "NM_000546.6")
        self.assertEqual(messages, [SearchMessage('Using transcript "NM_000546" version "6" instead of provided version "12"',
                                                  severity=LogLevel.ERROR, substituted=True)])

    def test_gene_symbol_message_plurals(self):
        single = _gene_symbol_search_message({"P53 (HGNC)"}, {"TP53"})
        self.assertIn("via alias: P53 (HGNC).", single)
        self.assertIn("Gene symbol TP53 has non-MANE", single)

        multiple = _gene_symbol_search_message({"A (HGNC)", "B (HGNC)"}, {"X", "Y"})
        self.assertIn("via aliases: ", multiple)
        self.assertRegex(multiple, r"Gene symbols (X, Y|Y, X) have non-MANE")

    def test_gene_symbol_message_mane_only(self):
        self.assertEqual(_gene_symbol_search_message(set(), set()),
                         "Searching without transcript is only supported for coordinates resolved against MANE transcripts.")

    def test_nomenclature_error_message(self):
        self.assertEqual(_hgvs_nomenclature_error_message("NM_000546.6:c.215C>G", "Readable problem"), "Readable problem")
        self.assertEqual(_hgvs_nomenclature_error_message("NM_000546.6:c.215C>Gxx", "unexpected char 'x'"),
                         'Invalid HGVS cDNA allele : "c.215C>Gxx"')


class ContigNotInBuildMessageTestCase(TestCase):

    def test_unknown_contig(self):
        message = _contig_not_in_build_message("NC_999999.1")
        self.assertEqual(message.message, "contig_accession='NC_999999.1' not in enabled builds: "
                                          f"{GenomeBuild.get_enabled_builds_comma_separated_string()}")
        self.assertEqual(message.severity, LogLevel.WARNING)

    @override_settings(ANNOTATION={"GRCh37": {"enabled": True}})
    def test_contig_in_searched_build_is_silent(self):
        self.assertIsNone(_contig_not_in_build_message("NC_000017.10"))

    @override_settings(ANNOTATION={"GRCh37": {"enabled": True}})
    def test_contig_in_unannotated_build(self):
        message = _contig_not_in_build_message("NC_000017.11")
        self.assertEqual(message.message, "contig_accession='NC_000017.11' is from builds GRCh38 which are not enabled/annotated.")


class SearchHGVSUsingGeneSymbolTestCase(TestCase):
    """ MANE transcripts resolving to different variants must each keep only their own messages (#1950) """

    def test_messages_stay_with_their_result(self):
        grch38 = GenomeBuild.grch38()

        def fake_search_hgvs(transcript_hgvs, *args):
            accession = transcript_hgvs.split(":")[0]
            identifier = "variant_b" if accession == "NM_B.1" else "variant_a"
            preview = PreviewData(category="Variant", identifier=identifier)
            yield SearchResult(preview=preview, messages=[SearchMessage(f"resolved {accession}")])

        candidates = {"ENST_A.1": "MANE Select", "NM_A.1": "MANE Select", "NM_B.1": "MANE Plus Clinical"}
        with mock.patch("snpdb.signals.variant_search._search_hgvs", side_effect=fake_search_hgvs):
            results = list(_search_hgvs_using_gene_symbol(candidates, "c.1A>G", None, grch38, Variant.objects.none()))

        messages_by_identifier = {r.preview.identifier: [m.message for m in r.messages] for r in results}
        self.assertEqual(messages_by_identifier, {
            "variant_a": ["Results for: ENST_A.1 (MANE Select), NM_A.1 (MANE Select)",
                          "resolved ENST_A.1", "resolved NM_A.1"],
            "variant_b": ["Results for: NM_B.1 (MANE Plus Clinical)", "resolved NM_B.1"],
        })
