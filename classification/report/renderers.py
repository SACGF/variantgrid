"""
One HTML document, and the PDF, Word file and JSON derived from it.

The lab's whole report design is the one Django template, rendered server side, so the preview the
scientist reads, the PDF filed with the case and the Word file they edit cannot disagree. The JSON
comes off the same context, so the structured record and the human document describe the same
classifications at the same versions - built by whichever app owns that template's shape, else as
the canonical context dump.

Entry point: render_case_report(template, context, title) -> RenderedCaseReport.
"""
import json
import re
from dataclasses import dataclass
from typing import Optional

from django.core.serializers.json import DjangoJSONEncoder
from html2docx import html2docx

from classification.models.classification_report_models import (
    ClassificationReportTemplate,
    case_report_json_signal,
)
from classification.report.template_validation import html_to_pdf, render_html


@dataclass
class RenderedCaseReport:
    html: str
    pdf: bytes
    docx: bytes
    json_output: dict
    # The context as stored - JSON safe, so the same documents can be re-rendered after a template fix
    context_snapshot: dict


def json_safe(context: dict) -> dict:
    """ Dates and decimals the way a JSONField will take them back """
    return json.loads(json.dumps(context, cls=DjangoJSONEncoder))


def render_json(report_template: ClassificationReportTemplate, context: dict) -> dict:
    """ The report's structured record: the app that owns this template's shape, else the canonical
        context dump - so there is always one.

        A JSON another system parses is an interface, and keeping it in step with that system is
        code rather than config a lab edits, so it is written in Python by the app that has to.
        Anything raised there fails the build rather than quietly falling through to a differently
        shaped document. """
    for _receiver, result in case_report_json_signal.send(sender=ClassificationReportTemplate,
                                                          report_template=report_template,
                                                          context=context):
        if result is not None:
            return result
    return json_safe(context)


# html2docx has no notion of a document head: it prints the contents of <style> and <script> as the
# first paragraph of the Word file. They are the print design, not the report, so they come out first
NON_CONTENT_PATTERN = re.compile(r"<(head|style|script)\b.*?</\1\s*>", re.DOTALL | re.IGNORECASE)


def render_docx(html: str, title: str) -> bytes:
    return html2docx(NON_CONTENT_PATTERN.sub("", html), title=title).getvalue()


def render_case_report(template: ClassificationReportTemplate, context: dict,
                       title: Optional[str] = None) -> RenderedCaseReport:
    html = render_html(template.case_template, context)
    return RenderedCaseReport(
        html=html,
        pdf=html_to_pdf(html),
        docx=render_docx(html, title or template.name),
        json_output=render_json(template, context),
        context_snapshot=json_safe(context),
    )
