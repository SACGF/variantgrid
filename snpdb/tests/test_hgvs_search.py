from unittest import mock

from django.test import SimpleTestCase, TestCase, override_settings

from library.enums.log_level import LogLevel
from library.preview_request import PreviewData
from snpdb.models import GenomeBuild, Variant
from snpdb.search import SearchMessage, SearchResult
from snpdb.signals.hgvs_search import (
    _contig_not_in_build_message,
    _gene_symbol_search_message,
    _hgvs_nomenclature_error_message,
    _search_hgvs_using_gene_symbol,
    _transcript_substitution_messages,
)


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
        with mock.patch("snpdb.signals.hgvs_search._search_hgvs", side_effect=fake_search_hgvs):
            results = list(_search_hgvs_using_gene_symbol(candidates, "c.1A>G", None, grch38, Variant.objects.none()))

        messages_by_identifier = {r.preview.identifier: [m.message for m in r.messages] for r in results}
        self.assertEqual(messages_by_identifier, {
            "variant_a": ["Results for: ENST_A.1 (MANE Select), NM_A.1 (MANE Select)",
                          "resolved ENST_A.1", "resolved NM_A.1"],
            "variant_b": ["Results for: NM_B.1 (MANE Plus Clinical)", "resolved NM_B.1"],
        })
