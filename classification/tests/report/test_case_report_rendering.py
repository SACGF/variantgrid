import json

from django.core.exceptions import ValidationError
from django.test import TestCase

from classification.enums.classification_enums import SomaticClinicalSignificance
from classification.models import ClassificationReportTemplate
from classification.report.case_report_context import (
    ReportContext,
    build_gene_groups,
    build_kind_groups,
    build_tier_groups,
    context_as_dict,
    sort_report_variants,
)
from classification.report.default_templates import generic_case_template, generic_json_template
from classification.report.renderers import render_case_report
from classification.tests.report.fake_report_variants import fake_report_variant

TIER_1 = SomaticClinicalSignificance.TIER_1

JSON_TEMPLATE = """{% load js_tags %}
{"Variants": [{% for v in variants %}
  {"Gene": {{ v.gene_label|jsonify }},
   "c": {{ v.evidence.c_hgvs.value|jsonify }},
   "AMPTier": "Tier {{ v.amp_tier }}",
   "Report": "{% if v.reported %}Y{% else %}N{% endif %}"}{% if not forloop.last %},{% endif %}
{% endfor %}]}"""


def _report_context(variants) -> dict:
    variants = sort_report_variants(variants)
    return context_as_dict(ReportContext(source_level="S", variants=variants,
                                         kind_groups=build_kind_groups(variants),
                                         tier_groups=build_tier_groups(variants),
                                         gene_groups=build_gene_groups(variants)))


class CaseReportRenderingTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.template = ClassificationReportTemplate.objects.create(
            name="test case template",
            case_template=generic_case_template(),
            json_template=JSON_TEMPLATE)
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
        """ "Report: N" means seen but not printed - it stays in the structured record so the
            downstream system knows the case was looked at """
        rendered = render_case_report(self.template,
                                      _report_context([self.reported, self.unreported]))

        self.assertNotIn("c.222C&gt;T", rendered.html)
        by_gene = {v["Gene"]: v for v in rendered.json_output["Variants"]}
        self.assertEqual(by_gene["AGENE"]["Report"], "Y")
        self.assertEqual(by_gene["BGENE"]["Report"], "N")
        self.assertEqual(by_gene["BGENE"]["c"], "c.222C>T")


class JsonTemplateTest(TestCase):

    def test_a_blank_json_template_dumps_the_context(self):
        template = ClassificationReportTemplate.objects.create(
            name="blank json", case_template=generic_case_template(),
            json_template=generic_json_template())
        context = _report_context([fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"],
                                                       c_hgvs="c.111A>G")])

        rendered = render_case_report(template, context)

        self.assertEqual(rendered.json_output["variants"][0]["gene_symbol"], "AGENE")
        self.assertEqual(rendered.json_output, rendered.context_snapshot)

    def test_a_template_that_renders_invalid_json_is_rejected(self):
        template = ClassificationReportTemplate(name="broken json",
                                                json_template='{"Variants": [,]}')
        with self.assertRaises(ValidationError) as raised:
            template.full_clean()
        self.assertIn("json_template", raised.exception.message_dict)

    def test_a_valid_json_template_saves(self):
        template = ClassificationReportTemplate(name="fine json", json_template=JSON_TEMPLATE)
        template.full_clean()  # raises if not

    def test_a_case_template_that_cannot_be_a_pdf_is_rejected(self):
        template = ClassificationReportTemplate(name="broken case",
                                                case_template="{% for v in variants %}")
        with self.assertRaises(ValidationError) as raised:
            template.full_clean()
        self.assertIn("case_template", raised.exception.message_dict)

    def test_the_shipped_default_templates_are_valid(self):
        template = ClassificationReportTemplate(name="shipped default",
                                                case_template=generic_case_template(),
                                                json_template=generic_json_template())
        template.full_clean()  # raises if not


class ContextSnapshotTest(TestCase):

    def test_the_snapshot_is_json_serialisable(self):
        """ context_snapshot is a JSONField, and the live context carries dates and datetimes """
        context = _report_context([fake_report_variant("AGENE", tier=TIER_1, amp_levels=["A"])])
        template = ClassificationReportTemplate.objects.create(
            name="snapshot", case_template="<html><body>x</body></html>")

        rendered = render_case_report(template, context)

        json.dumps(rendered.context_snapshot)  # raises if not
