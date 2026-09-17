import re
from typing import Optional

from bs4 import BeautifulSoup
from django.contrib.auth.models import User
from django.http import HttpResponse
from django.template import engines
from django.utils import timezone

from classification.models import ClassificationReportTemplate
from classification.models.classification import ClassificationModification
from classification.report.case_report_context import (
    build_gene_groups,
    build_report_variants,
    evidence_row_data,
)

UNSUBMITTED_CHANGES_WARNING = (
    '<div style="margin: 8px; padding: 8px 12px; border: 1px solid #f5c6cb; border-radius: 4px;'
    ' background-color: #f8d7da; color: #721c24; font-family: sans-serif; font-size: 14px;">{message}</div>'
)


class ClassificationReport:
    """
    Formats using report for the corresponding lab.

    One record fills "record" the way it always has. A multi-variant report (the sample/patient page's
    Classify & Report tab) passes the rest in as well - "classifications" is every record in report order and
    "gene_groups" the same records grouped by gene symbol - so single-variant templates keep working unchanged.
    @see classification/report/case_report_context.py, which builds the rows for both this and the case report.
    """

    def __init__(self, classification: ClassificationModification, user: User,
                 classifications: Optional[list[ClassificationModification]] = None,
                 report_template: Optional[ClassificationReportTemplate] = None,
                 extra_context: Optional[dict] = None):
        self.classification = classification
        self.classifications = classifications if classifications is not None else [classification]
        self.report_template = report_template
        self.extra_context = extra_context or {}  # eg the case header a multi-variant report opens with
        self.user = user

    def get_template(self):
        report = self.report_template or ClassificationReportTemplate.preferred_template_for(self.classification)
        template_str = report.template or 'No report template has been configured'
        django_engine = engines['django']
        return django_engine.from_string(template_str)

    def context(self) -> dict:
        """ The single record report and the case report build their rows from the same ReportVariant,
            so "this variant's row" means one thing - `record` / `classifications` / `gene_groups` are
            the shape the existing templates read """
        variants = build_report_variants(self.classifications, self.user)
        by_pk = {variant.modification.pk: variant for variant in variants}

        record = by_pk[self.classification.pk].evidence if self.classification.pk in by_pk \
            else self.row_data(self.classification)
        return {
            **self.extra_context,
            'record': record,
            'classifications': [variant.evidence for variant in variants],
            'gene_groups': build_gene_groups(variants),
        }

    def _unsubmitted_changes_warning(self) -> Optional[str]:
        """ Reports are rendered from submitted versions, so say so when the form is showing newer data """
        stale = [record for record in self.classifications if not record.is_last_edited]
        if not stale:
            return None

        if len(self.classifications) == 1:
            submitted = timezone.localtime(stale[0].created).strftime("%d %b %Y %H:%M")
            message = (f"This report was generated from the version submitted on {submitted}. "
                       "Changes made since then have not been submitted and do not appear below.")
        else:
            message = (f"{len(stale)} of {len(self.classifications)} classifications have unsubmitted changes - "
                       "this report was generated from their submitted versions.")
        return UNSUBMITTED_CHANGES_WARNING.format(message=message)

    @staticmethod
    def _insert_warning(content: str, warning: str) -> str:
        """ The templates are lab maintained config, so the warning goes in the rendered page rather than in them """
        if body := re.search(r"<body[^>]*>", content, re.IGNORECASE):
            return content[:body.end()] + warning + content[body.end():]
        return warning + content

    def serve(self):
        template = self.get_template()
        content = template.render(self.context())
        if warning := self._unsubmitted_changes_warning():
            content = self._insert_warning(content, warning)
        response = HttpResponse(content=content, content_type='text/html')
        return response

    def row_data(self, record: ClassificationModification) -> dict:
        return evidence_row_data(record, self.user)

    def get_unknown_evidence(self) -> list[str]:
        row_data = self.row_data(self.classification)
        template = self.get_template()
        content = template.render({'record': row_data})
        return self._get_unknown_evidence(content, row_data)

    @staticmethod
    def _get_unknown_evidence(content, row_data) -> list[str]:
        soup = BeautifulSoup(content, "html.parser")
        unknown_evidence = set()
        for tag in soup.find_all(attrs={":evidence": True}):
            evidence = tag[":evidence"]
            if re.match(r"[{(.]", evidence):
                continue
            if evidence not in row_data:
                unknown_evidence.add(evidence)
        return sorted(unknown_evidence)
