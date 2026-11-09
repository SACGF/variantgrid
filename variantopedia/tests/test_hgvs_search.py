from unittest import mock

from django.conf import settings
from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase, override_settings

from annotation.fake_data import get_fake_annotation_version
from genes.fake_data import create_gata2_transcript_version
from genes.models import MissingTranscript
from library.enums.log_level import LogLevel
from library.preview_request import PreviewData
from snpdb.models import GenomeBuild, Variant
from snpdb.search import SearchMessage, SearchResult, search_data
from variantopedia.signals.hgvs_search import (
    _contig_not_in_build_message,
    _gene_symbol_search_message,
    _hgvs_nomenclature_error_message,
    _hgvs_strings_for_variant_label,
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
        with mock.patch("variantopedia.signals.hgvs_search._search_hgvs", side_effect=fake_search_hgvs):
            results = list(_search_hgvs_using_gene_symbol(candidates, "c.1A>G", None, grch38, Variant.objects.none()))

        messages_by_identifier = {r.preview.identifier: [m.message for m in r.messages] for r in results}
        self.assertEqual(messages_by_identifier, {
            "variant_a": ["Results for: ENST_A.1 (MANE Select), NM_A.1 (MANE Select)",
                          "resolved ENST_A.1", "resolved NM_A.1"],
            "variant_b": ["Results for: NM_B.1 (MANE Plus Clinical)", "resolved NM_B.1"],
        })

    def test_transcript_missing_for_every_candidate(self):
        """ None of the symbol's transcripts being in the build isn't a resolution error to report (#2692) """
        def missing_search_hgvs(transcript_hgvs, *args):
            raise MissingTranscript(f"{transcript_hgvs.split(':')[0]} is not available for T2T-CHM13v2.0")

        with mock.patch("variantopedia.signals.hgvs_search._search_hgvs", side_effect=missing_search_hgvs):
            with self.assertRaises(MissingTranscript):
                list(_search_hgvs_using_gene_symbol({"NM_A.1": "MANE Select"}, "c.1A>G", None,
                                                    GenomeBuild.grch38(), Variant.objects.none()))


@override_settings(PREFER_ALLELE_LINKS=False)
class SearchHGVSTranscriptMissingFromBuildTestCase(TestCase):
    """ A transcript in one build only (T2T has no ClinGen to fall back on) is reported against the build it is
        missing from, not sent to Rollbar (#2692) """
    HGVS = "NM_001145661.2(GATA2):c.1121G>A"
    MISSING_MESSAGE = "NM_001145661.2 is not available for T2T-CHM13v2.0"

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="testuser")[0]
        cls.grch37 = GenomeBuild.grch37()
        cls.t2t = GenomeBuild.t2tv2()
        get_fake_annotation_version(cls.grch37)
        get_fake_annotation_version(cls.t2t)
        create_gata2_transcript_version(cls.grch37)

    def _search(self, *annotated_builds: GenomeBuild):
        annotated_names = {genome_build.name for genome_build in annotated_builds}
        annotation = {name: {**build_settings, "enabled": name in annotated_names}
                      for name, build_settings in settings.ANNOTATION.items()}
        with self.settings(ANNOTATION=annotation), \
                mock.patch("variantopedia.signals.hgvs_search.report_exc_info") as report_exc_info:
            search_results = search_data(self.user, self.HGVS, False)
        report_exc_info.assert_not_called()
        return search_results

    def test_info_when_another_build_resolved(self):
        search_results = self._search(self.grch37, self.t2t)
        self.assertFalse([r for r in search_results.results if r.preview.is_error])
        missing = [m for m in search_results.messages_overall if m.message == self.MISSING_MESSAGE]
        self.assertEqual([(m.severity, m.genome_builds) for m in missing], [(LogLevel.INFO, [self.t2t])])

    def test_error_when_no_build_resolved(self):
        search_results = self._search(self.t2t)
        errors = [(m.message, r.preview.genome_builds) for r in search_results.results for m in r.messages]
        self.assertEqual(errors, [(self.MISSING_MESSAGE, {self.t2t})])


class VariantLabelHGVSTestCase(TestCase):
    """ The Variant column's label leads HGVS with a gene, the genes an SV spans or the contig (#2696) """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        create_gata2_transcript_version(cls.grch37)

    def _hgvs(self, label: str) -> list[str]:
        return _hgvs_strings_for_variant_label(label, self.grch37)

    def test_genomic_change_goes_on_the_contig(self):
        self.assertEqual(self._hgvs("chr3 g.128198980A>G"), ["NC_000003.11:g.128198980A>G"])
        self.assertEqual(self._hgvs("chrMT m.5782_13922del"), ["NC_012920.1:m.5782_13922del"])
        self.assertEqual(self._hgvs("chr3 +5 genes g.128198980_128199980del"), ["NC_000003.11:g.128198980_128199980del"])

    def test_genomic_change_goes_on_the_genes_contig(self):
        self.assertEqual(self._hgvs("GATA2, RPN1 +3 genes g.128198980_128199980del"),
                         ["NC_000003.11:g.128198980_128199980del"])

    def test_transcript_change_keeps_the_first_gene(self):
        self.assertEqual(self._hgvs("GATA2, RPN1 c.1113dup"), ["GATA2:c.1113dup"])

    def test_hgvs_that_is_not_a_label_is_unchanged(self):
        for hgvs in ["GATA2 c.1113dup", "NM_001145661.2(GATA2):c.1113dup", "NM_001145661.2:c.1113dup"]:
            self.assertEqual(self._hgvs(hgvs), [hgvs])
