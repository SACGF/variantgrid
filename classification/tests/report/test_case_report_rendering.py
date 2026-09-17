import io
import json
import zipfile

from django.core.exceptions import ValidationError
from django.test import TestCase

from classification.enums.classification_enums import SomaticClinicalSignificance
from classification.models import ClassificationReportTemplate
from classification.models.classification_report_models import case_report_json_signal
from classification.report.case_report_context import (
    ReportContext,
    build_gene_groups,
    build_kind_groups,
    build_tier_groups,
    context_as_dict,
    sort_report_variants,
)
from classification.report.default_templates import generic_case_template
from classification.report.renderers import (
    DRAFT_BANNER_TEXT,
    PRE_FONT_NAME,
    render_case_report,
    render_docx,
)
from classification.tests.report.fake_report_variants import fake_report_variant

TIER_1 = SomaticClinicalSignificance.TIER_1


def _report_context(variants, draft: bool = True) -> dict:
    variants = sort_report_variants(variants)
    return context_as_dict(ReportContext(source_level="S", variants=variants,
                                         kind_groups=build_kind_groups(variants),
                                         tier_groups=build_tier_groups(variants),
                                         gene_groups=build_gene_groups(variants),
                                         draft=draft))


# The class the renderer strips for the Word file - the CSS rule for it is in the head either way
WATERMARK_ELEMENT = '<div class="draft-watermark">'


def _docx_xml(docx: bytes) -> str:
    """ What Word will actually read - the only way to see what html2docx made of the HTML """
    with zipfile.ZipFile(io.BytesIO(docx)) as archive:
        return archive.read("word/document.xml").decode("utf-8")


class CaseReportRenderingTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.template = ClassificationReportTemplate.objects.create(
            name="test case template",
            case_template=generic_case_template())
        cls.reported = fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], vaf=0.42,
                                           c_hgvs="c.111A>G", pk=1)
        cls.unreported = fake_report_variant("BGENE", tier=TIER_1, amp_levels=["A"], vaf=0.2,
                                             c_hgvs="c.222C>T", reported=False, pk=2)

    def test_the_documents_all_come_from_the_one_html(self):
        rendered = render_case_report(self.template,
                                      _report_context([self.reported, self.unreported]))

        self.assertEqual(rendered.html.count("c.111A&gt;G"), 1)
        self.assertTrue(rendered.pdf.startswith(b"%PDF"))
        self.assertTrue(rendered.docx)

    def test_an_unreported_variant_is_off_the_document_and_on_the_record(self):
        """ Seen but not printed - it stays in the structured record so the downstream system knows
            the case was looked at """
        rendered = render_case_report(self.template,
                                      _report_context([self.reported, self.unreported]))

        self.assertNotIn("c.222C&gt;T", rendered.html)
        by_gene = {v["gene_symbol"]: v for v in rendered.json_output["variants"]}
        self.assertTrue(by_gene["AGENE"]["reported"])
        self.assertFalse(by_gene["BGENE"]["reported"])
        self.assertEqual(by_gene["BGENE"]["evidence"]["c_hgvs"]["value"], "c.222C>T")


class ReportJsonTest(TestCase):

    def _context(self) -> dict:
        return _report_context([fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"],
                                                    c_hgvs="c.111A>G")])

    def test_with_no_app_answering_the_context_is_the_record(self):
        template = ClassificationReportTemplate.objects.create(
            name="no app json", case_template=generic_case_template())

        rendered = render_case_report(template, self._context())

        self.assertEqual(rendered.json_output["variants"][0]["gene_symbol"], "AGENE")
        self.assertEqual(rendered.json_output, rendered.context_snapshot)

    def test_an_app_that_owns_the_shape_writes_the_json(self):
        """ A JSON another system parses is an interface, so the app keeping it in step with that
            system writes it - SA Path's TSO 500 JSON arrives this way """
        template = ClassificationReportTemplate.objects.create(
            name="app owned json", case_template=generic_case_template())

        def json_receiver(sender, report_template, context, **kwargs):
            return {"built": "in python", "for": report_template.pk}

        case_report_json_signal.connect(json_receiver)
        try:
            rendered = render_case_report(template, self._context())
        finally:
            case_report_json_signal.disconnect(json_receiver)

        self.assertEqual(rendered.json_output, {"built": "in python", "for": "app owned json"})

    def test_a_case_template_that_cannot_be_a_pdf_is_rejected(self):
        template = ClassificationReportTemplate(name="broken case",
                                                case_template="{% for v in variants %}")
        with self.assertRaises(ValidationError) as raised:
            template.full_clean()
        self.assertIn("case_template", raised.exception.message_dict)

    def test_the_shipped_default_template_is_valid(self):
        template = ClassificationReportTemplate(name="shipped default",
                                                case_template=generic_case_template())
        template.full_clean()  # raises if not


class ContextSnapshotTest(TestCase):

    def test_the_snapshot_is_json_serialisable(self):
        """ context_snapshot is a JSONField, and the live context carries dates and datetimes """
        context = _report_context([fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"])])
        template = ClassificationReportTemplate.objects.create(
            name="snapshot", case_template="<html><body>x</body></html>")

        rendered = render_case_report(template, context)

        json.dumps(rendered.context_snapshot)  # raises if not


class DraftWatermarkTest(TestCase):
    """ A draft has to be unmistakable on paper, and the two document formats say so differently:
        the PDF gets the template's watermark element, the Word file a banner the renderer adds
        because html2docx has nothing that can be positioned or rotated """

    @classmethod
    def setUpTestData(cls):
        cls.template = ClassificationReportTemplate.objects.create(
            name="draft watermark", case_template=generic_case_template())
        cls.variants = [fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"], vaf=0.42)]

    def test_a_draft_is_watermarked_in_the_pdf_and_bannered_in_the_word_file(self):
        rendered = render_case_report(self.template, _report_context(self.variants, draft=True))

        self.assertIn(WATERMARK_ELEMENT, rendered.html)
        docx_xml = _docx_xml(rendered.docx)
        self.assertIn(DRAFT_BANNER_TEXT, docx_xml)
        # The watermark element itself is print design - it would arrive as a stray DRAFT paragraph
        self.assertEqual(docx_xml.count("DRAFT"), 1)

    def test_the_issued_report_carries_neither(self):
        rendered = render_case_report(self.template, _report_context(self.variants, draft=False))

        self.assertNotIn(WATERMARK_ELEMENT, rendered.html)
        self.assertNotIn("DRAFT", _docx_xml(rendered.docx))


class ResultsSummaryDocxTest(TestCase):
    """ A Results Summary is column aligned text, and html2docx gives a <pre> block its spaces and
        line breaks but no font at all - in Word's proportional default nothing lines up """

    PRE_HTML = ('<html><body><pre class="results-summary">Somatic Variants   Tier   Gene\n'
                '                   IA     BRAF\n</pre></body></html>')

    def test_a_pre_block_reaches_word_monospaced_and_line_for_line(self):
        docx_xml = _docx_xml(render_docx(self.PRE_HTML, title="pre"))

        self.assertIn(f'w:ascii="{PRE_FONT_NAME}"', docx_xml)
        self.assertIn("<w:br/>", docx_xml)
        self.assertIn('<w:t xml:space="preserve">                   IA     BRAF</w:t>', docx_xml)
