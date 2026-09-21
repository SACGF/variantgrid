"""
A built case report's own actions: the four downloads, Finalise, New version, Rebuild documents and
the LIS fields (external_report_id / report_date).

These take a CaseReport rather than a case, which is why they live here rather than on the analysis
app's Classify & Report tab - the report already knows which case it is about, and
CaseReport.can_view / can_write decide who may act. Building one needs the case's classifications,
so that stays in analysis/views/views_classify_report.py.

Every document is served through a view rather than a media URL: MEDIA_ROOT is not permission
checked, and a case report is a patient's.
"""
import json

from django.core.exceptions import PermissionDenied
from django.http import FileResponse, Http404, HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404
from django.utils.dateparse import parse_date
from django.views.decorators.clickjacking import xframe_options_sameorigin
from django.views.decorators.http import require_POST

from classification.models import CaseReport
from classification.report.case_report_builder import (
    FinaliseResult,
    build_new_version,
    finalise_case_report,
    rebuild_documents,
    stamp_report_onto_classifications,
)

HTML_FORMAT = "html"
JSON_FORMAT = "json"
FILE_FORMATS = {
    "pdf": "pdf_file",
    "docx": "docx_file",
}


def _get_case_report(user, case_report_id: int, write: bool = False) -> CaseReport:
    case_report = get_object_or_404(CaseReport, pk=case_report_id)
    if write:
        case_report.check_can_write(user)
    else:
        case_report.check_can_view(user)
    return case_report


def _attachment(content: bytes, content_type: str, filename: str) -> HttpResponse:
    response = HttpResponse(content=content, content_type=content_type)
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


def case_report_download(request, case_report_id: int, document_format: str):
    """ The report as the scientist filed it - checked against the report's lab and the case, since
        MEDIA_ROOT has no permissions of its own """
    case_report = _get_case_report(request.user, case_report_id)
    stem = f"case_report_{case_report.pk}"

    if document_format == HTML_FORMAT:
        return _attachment(case_report.html.encode("utf-8"), "text/html", f"{stem}.html")
    if document_format == JSON_FORMAT:
        content = json.dumps(case_report.json_output, indent=2).encode("utf-8")
        return _attachment(content, "application/json", f"{stem}.json")

    field_name = FILE_FORMATS.get(document_format)
    if field_name is None:
        raise Http404(f"Unknown report format '{document_format}'")
    file_field = getattr(case_report, field_name)
    if not file_field:
        raise Http404(f"Case report {case_report.pk} has no {document_format} document")
    return FileResponse(file_field.open("rb"), as_attachment=True, filename=f"{stem}.{document_format}")


@xframe_options_sameorigin
def view_case_report(request, case_report_id: int) -> HttpResponse:
    """ The stored HTML, as the preview the scientist reads - framed by the built-report modal, which the
        site-wide X-Frame-Options DENY would otherwise blank """
    case_report = _get_case_report(request.user, case_report_id)
    return HttpResponse(content=case_report.html)


def _finalise_response(result: FinaliseResult) -> JsonResponse:
    return JsonResponse({
        "case_report_id": result.case_report.pk,
        "status": result.case_report.get_status_display(),
        "stamped": [c.cr_lab_id for c in result.stamped],
        # Named rather than published - the report never pushes out edits nobody has submitted
        "unsubmitted": [c.cr_lab_id for c in result.unsubmitted],
    })


@require_POST
def case_report_finalise(request, case_report_id: int) -> JsonResponse:
    case_report = _get_case_report(request.user, case_report_id, write=True)
    return _finalise_response(finalise_case_report(case_report, request.user))


@require_POST
def case_report_rebuild(request, case_report_id: int) -> JsonResponse:
    """ Re-render the documents from the context the report was built over - a template fix """
    case_report = _get_case_report(request.user, case_report_id, write=True)
    if not case_report.is_editable:
        raise PermissionDenied(f"Case report {case_report.pk} is {case_report.get_status_display()} - "
                               "the documents that went out with the case cannot be re-rendered")
    rebuild_documents(case_report)
    return JsonResponse({"case_report_id": case_report.pk})


@require_POST
def case_report_new_version(request, case_report_id: int) -> JsonResponse:
    case_report = _get_case_report(request.user, case_report_id, write=True)
    new_version = build_new_version(case_report, request.user)
    return JsonResponse({"case_report_id": new_version.pk, "supersedes": case_report.pk})


@require_POST
def case_report_lis_details(request, case_report_id: int) -> JsonResponse:
    """ The LIS's report ID and date, entered once it has issued them. On a report that is already
        final this stamps report_id onto the pinned classifications without making a new version """
    case_report = _get_case_report(request.user, case_report_id, write=True)
    case_report.external_report_id = request.POST.get("external_report_id") or None
    case_report.report_date = parse_date(request.POST.get("report_date") or "")
    case_report.save()
    result = stamp_report_onto_classifications(case_report, request.user)
    return _finalise_response(result)
