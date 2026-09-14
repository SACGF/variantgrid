"""
Building, rebuilding and finalising a CaseReport.

A build renders the context once and keeps all of it: classification/report/renderers.py turns the
one HTML into the PDF and the Word file, and `context_snapshot` is what `Rebuild documents`
re-renders from after a template fix, so a report's numbers never move under it. `New version`
instead builds a fresh context from each pinned classification's current published version and
supersedes the old report.

Finalising is the other half. The report becomes the copy that went out with the case, and the three
keys that say so - report_date, variant_reported and report_id - are stamped onto each pinned
classification and re-published. Going FINAL sends case_report_finalised_signal, which is how a
deployment specific app files the report somewhere else. Every stamp is idempotent (patch_value
writes nothing when the values already match), and a classification with unsubmitted edits is left
alone and named, because publishing it would push someone's work in progress out under the report's
name.

Entry points: build_case_report, preview_case_report_html, rebuild_documents, build_new_version,
finalise_case_report, stamp_report_onto_classifications.
"""
from dataclasses import dataclass, field
from typing import Optional

from django.contrib.auth.models import User
from django.core.files.base import ContentFile
from django.utils import timezone

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_report_models import (
    CaseReport,
    CaseReportClassification,
    CaseReportStatus,
    ClassificationReportTemplate,
    case_report_finalised_signal,
)
from classification.report.case_report_context import (
    NOT_REPORTED,
    build_report_context,
    case_report_as_dict,
    context_as_dict,
)
from classification.report.renderers import render_case_report
from classification.report.template_validation import render_html
from snpdb.models import Lab

# What variant_reported becomes for a variant the report prints that has never said which kind of
# finding it is. A record that already says keeps what it says - the report only ever settles
# "on the document" against "seen, not reported"
DEFAULT_REPORTED = "primary_finding"

PDF_FILENAME = "case_report.pdf"
DOCX_FILENAME = "case_report.docx"


@dataclass
class FinaliseResult:
    """ What finalising did. `unsubmitted` are the records left alone - the finalise dialog lists
        them, since publishing one would push out edits nobody has submitted """
    case_report: CaseReport
    stamped: list[Classification] = field(default_factory=list)
    unsubmitted: list[Classification] = field(default_factory=list)


def _pinned_rows(case_report: CaseReport) -> list[CaseReportClassification]:
    """ Every pinned row with its classification loaded - the builder reads one off each """
    return list(case_report.casereportclassification_set
                .select_related("classification_modification__classification"))


def _save_document(case_report: CaseReport, field_name: str, filename: str, content: bytes):
    """ Rebuilding replaces the file rather than leaving a suffixed copy beside it """
    file_field = getattr(case_report, field_name)
    if file_field:
        file_field.delete(save=False)
    file_field.save(filename, ContentFile(content), save=False)


def _render_and_save(case_report: CaseReport, context: dict) -> CaseReport:
    rendered = render_case_report(case_report.template, context, title=str(case_report.source))
    case_report.html = rendered.html
    case_report.json_output = rendered.json_output
    case_report.context_snapshot = rendered.context_snapshot
    _save_document(case_report, "pdf_file", PDF_FILENAME, rendered.pdf)
    _save_document(case_report, "docx_file", DOCX_FILENAME, rendered.docx)
    case_report.save()
    return case_report


def build_case_report(user: User, template: ClassificationReportTemplate, lab: Lab,
                      source_level: str, source,
                      modifications: list[ClassificationModification],
                      reported_by_pk: Optional[dict[int, bool]] = None,
                      summary: str = "", case_values: Optional[dict] = None,
                      supersedes: Optional[CaseReport] = None) -> CaseReport:
    """ One run of the template over the case: the CaseReport, its pinned rows in report order, and
        the four documents. The report is created first because its pk is in the media path """
    case_report = CaseReport.objects.create(template=template, lab=lab, user=user,
                                            supersedes=supersedes, summary=summary,
                                            case_values=case_values or {},
                                            **CaseReport.source_kwargs(source_level, source))
    report_context = build_report_context(user, source_level, source, modifications, lab=lab,
                                          reported_by_pk=reported_by_pk, summary=summary,
                                          case_values=case_values, case_report=case_report)
    CaseReportClassification.objects.bulk_create([
        CaseReportClassification(case_report=case_report, classification_modification=variant.modification,
                                 order=order, reported=variant.reported)
        for order, variant in enumerate(report_context.variants)])
    _render_and_save(case_report, context_as_dict(report_context))

    if supersedes and supersedes.status != CaseReportStatus.SUPERSEDED:
        CaseReport.objects.filter(pk=supersedes.pk).update(status=CaseReportStatus.SUPERSEDED)
    return case_report


def preview_case_report_html(user: User, template: ClassificationReportTemplate,
                             source_level: str, source,
                             modifications: list[ClassificationModification],
                             reported_by_pk: Optional[dict[int, bool]] = None,
                             summary: str = "", case_values: Optional[dict] = None,
                             lab: Optional[Lab] = None) -> str:
    """ The build form's Preview - the document as it would read, with nothing saved """
    report_context = build_report_context(user, source_level, source, modifications, lab=lab,
                                          reported_by_pk=reported_by_pk, summary=summary,
                                          case_values=case_values)
    return render_html(template.case_template, context_as_dict(report_context))


def rebuild_documents(case_report: CaseReport) -> CaseReport:
    """ A template fix, re-rendered over the context the report was built from - the numbers a
        report has already quoted stay as they were. Only the report's own block is refreshed, so
        an external_report_id entered after the build reaches the documents """
    context = dict(case_report.context_snapshot or {})
    context["case_report"] = case_report_as_dict(case_report)
    return _render_and_save(case_report, context)


def _current_published(case_report: CaseReport) -> tuple[list[ClassificationModification], dict[int, bool]]:
    """ Each pinned classification at its current published version, keeping the Report Y/N the
        scientist chose. A record with nothing published stays at the version that was rendered """
    modifications = []
    reported_by_pk = {}
    for row in _pinned_rows(case_report):
        classification = row.classification_modification.classification
        modification = classification.last_published_version or row.classification_modification
        modifications.append(modification)
        reported_by_pk[modification.pk] = row.reported
    return modifications, reported_by_pk


def build_new_version(case_report: CaseReport, user: User) -> CaseReport:
    """ The same case, built again from what is published now - the old report is superseded rather
        than edited, so the copy that went out with the case is still there to read """
    modifications, reported_by_pk = _current_published(case_report)
    return build_case_report(user, case_report.template, case_report.lab, case_report.source_level,
                             case_report.source, modifications, reported_by_pk=reported_by_pk,
                             summary=case_report.summary, case_values=case_report.case_values,
                             supersedes=case_report)


def has_unsubmitted_edits(classification: Classification) -> bool:
    """ A version edited since it was last published is someone's work in progress """
    last_published = classification.last_published_version
    return last_published is None or not last_published.is_last_edited


def _variant_reported_value(classification: Classification, reported: bool) -> str:
    """ The key says which kind of finding it was, which is the curator's call - so a record that
        already names one keeps it, and the report only settles reported against not """
    if not reported:
        return NOT_REPORTED
    current = classification.get(SpecialEKeys.VARIANT_REPORTED)
    if current and current != NOT_REPORTED:
        return current
    return DEFAULT_REPORTED


def stamp_report_onto_classifications(case_report: CaseReport, user: User) -> FinaliseResult:
    """ report_date, variant_reported and report_id onto every pinned classification, then publish.
        patch_value writes nothing when the values already match, so re-finalising, rebuilding or
        re-entering the same LIS id creates no further modifications. The row then re-points at the
        new published version, which differs from the rendered one only in those keys, so the tab's
        stale check stays quiet and the pinned version is the one Shariant receives """
    result = FinaliseResult(case_report=case_report)
    for row in _pinned_rows(case_report):
        classification = row.classification_modification.classification
        if has_unsubmitted_edits(classification):
            result.unsubmitted.append(classification)
            continue

        patch = {SpecialEKeys.VARIANT_REPORTED: {"value": _variant_reported_value(classification, row.reported)}}
        if case_report.report_date:
            patch[SpecialEKeys.REPORT_DATE] = {"value": case_report.report_date.isoformat()}
        if case_report.external_report_id:
            patch[SpecialEKeys.REPORT_ID] = {"value": case_report.external_report_id}

        patch_response = classification.patch_value(patch, user=user,
                                                    source=SubmissionSource.VARIANT_GRID, save=True)
        if not patch_response.modified_keys:
            continue
        classification.publish_latest(user=user)
        if published := classification.last_published_version:
            CaseReportClassification.objects.filter(pk=row.pk).update(classification_modification=published)
        result.stamped.append(classification)
    return result


def finalise_case_report(case_report: CaseReport, user: User) -> FinaliseResult:
    """ The report becomes the copy that went out with the case, and says so on its classifications.
        Finalising an already final report re-stamps, which is how a LIS id entered later reaches
        the records without making a new version """
    newly_final = case_report.status == CaseReportStatus.DRAFT
    if newly_final:
        case_report.status = CaseReportStatus.FINAL
        if case_report.report_date is None:
            case_report.report_date = timezone.localdate()
        case_report.save()
    result = stamp_report_onto_classifications(case_report, user)
    if newly_final:
        # After the stamping, so a receiver filing the report elsewhere sees the finished state
        case_report_finalised_signal.send(sender=CaseReport, case_report=case_report, user=user)
    return result
