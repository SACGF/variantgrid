"""
Views for the classification grouping export ("Export NEW"): the settings page
(view_classification_grouping_export), the summary page that starts the download
(view_classification_grouping_export_redirect) and the export itself (serve_export).
"""
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render
from django.urls import reverse

from classification.lab_classification_stats import labs_with_classifications_qs
from classification.views.exports.classification_export_filter import InvalidExportParameter
from classification.views.exports.classification_export_formatter_vcf import VCFTargetSystem
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFilter, ClassificationGroupingExportFileSettings, ClassificationGroupingExportFormat
from classification.views.exports_grouping.classification_grouping_export_formatter_csv import \
    ClassificationGroupingExportFormatterCSV
from classification.views.exports_grouping.classification_grouping_export_formatter_json import \
    ClassificationGroupingExportFormatterJSON
from classification.views.exports_grouping.classification_grouping_export_formatter_vcf import \
    ClassificationGroupingExportFormatterVCF, VCFFormatDetails
from classification.views.exports_grouping.classification_grouping_export_process import \
    ClassificationGroupingExportProcess
from classification.views.exports_grouping.classification_grouping_export_formatter_franklin import \
    ClassificationGroupingExportFormatterFranklin, FranklinFormatDetails
from library.django_utils import get_url_from_view_path
from snpdb.enums import AlleleOriginBucket
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import GenomeBuild, Lab

FORMAT_LABELS = {
    "csv": "CSV",
    "vcf": "VCF",
    "json": "JSON",
    "franklin": "Franklin",
}

VCF_TARGET_SYSTEM_LABELS = {
    VCFTargetSystem.EMEDGENE: "Emedgene",
    VCFTargetSystem.VARSEQ: "VarSeq",
    VCFTargetSystem.GENERIC: "Generic",
}

LAB_MODE_LABELS = {
    "all": "All labs",
    "include-all-labs": "All labs",
    "exclude-my-labs": "Exclude my labs",
    "include-only-my-labs": "Include only my labs",
    "custom-labs": "Custom",
}


@dataclass(frozen=True)
class ClassificationGroupingExportSummary:
    format_label: str
    lab_selection: str
    filename: str
    query_string: str
    genome_build: Optional[GenomeBuild] = None
    vcf_target_system: Optional[str] = None
    included_labs: list[Lab] = field(default_factory=list)
    excluded_labs: list[Lab] = field(default_factory=list)
    allele_origin: Optional[AlleleOriginBucket] = None
    since: Optional[datetime] = None
    since_raw: Optional[str] = None
    rows_per_file: Optional[int] = None

    @property
    def since_is_days(self) -> bool:
        return bool(self.since_raw) and self.since_raw.isdigit()


def view_classification_grouping_export(request: HttpRequest) -> HttpResponseBase:
    user_labs = set(lab for lab in Lab.valid_labs_qs(request.user, admin_check=False))
    genome_builds = GenomeBuild.objects.all()
    return render(request, "classification/classification_grouping_export.html", {
        "genome_builds": genome_builds,
        "default_genome_build": GenomeBuildManager.get_current_genome_build(),
        "user_labs": user_labs,
        "sharing_labs": list(labs_with_classifications_qs(shared_only=True).filter(organization__active=True)
                             .select_related("organization").order_by("organization__name", "name")),
        "base_url": get_url_from_view_path(reverse('classification_grouping_export')),
        "redirect_url": get_url_from_view_path(reverse('classification_grouping_export_redirect')),
    })


def _export_format_from_request(request: HttpRequest,
                                export_filter: ClassificationGroupingExportFilter) -> ClassificationGroupingExportFormat:
    format_str = request.GET.get("format")
    match format_str:
        case "csv":
            return ClassificationGroupingExportFormatterCSV(classification_grouping_filter=export_filter)
        case "vcf":
            vcf_format_details = VCFFormatDetails(
                genome_build=GenomeBuild.get_name_or_alias(request.GET.get("genome_build")),
                target_system=VCFTargetSystem(request.GET.get("vcf_target_system"))
            )
            return ClassificationGroupingExportFormatterVCF(
                classification_grouping_filter=export_filter,
                vcf_format_details=vcf_format_details
            )
        case "json":
            return ClassificationGroupingExportFormatterJSON(classification_grouping_filter=export_filter)
        case "franklin":
            franklin_format_details = FranklinFormatDetails(
                genome_build=GenomeBuild.get_name_or_alias(request.GET.get("genome_build"))
            )
            return ClassificationGroupingExportFormatterFranklin(
                classification_grouping_filter=export_filter,
                franklin_formatter_details=franklin_format_details
            )
        case _:
            raise ValueError(f"Unsupported export format: {format_str}")


def _export_summary(request: HttpRequest) -> ClassificationGroupingExportSummary:
    export_filter = ClassificationGroupingExportFilter.from_request(request)
    export_settings = ClassificationGroupingExportFileSettings.from_request(request)
    export_format = _export_format_from_request(request, export_filter)

    genome_build = None
    vcf_target_system = None
    if isinstance(export_format, ClassificationGroupingExportFormatterVCF):
        genome_build = export_format.vcf_format_details.genome_build
        target_system = export_format.vcf_format_details.target_system
        vcf_target_system = VCF_TARGET_SYSTEM_LABELS.get(target_system, target_system.value)
    elif isinstance(export_format, ClassificationGroupingExportFormatterFranklin):
        genome_build = export_format.franklin_formatter_details.genome_build

    process = ClassificationGroupingExportProcess(
        classification_export_format=export_format,
        export_settings=export_settings,
        params=dict(request.GET),
        user=request.user
    )
    filename = process.filename(extension_override="zip" if export_settings.rows_per_file else None)

    def sorted_labs(labs) -> list[Lab]:
        return list(Lab.objects.filter(pk__in=[lab.pk for lab in labs]).select_related("organization")
                    .order_by("organization__name", "name"))

    return ClassificationGroupingExportSummary(
        format_label=FORMAT_LABELS[request.GET["format"]],
        lab_selection=LAB_MODE_LABELS[export_filter.lab_mode],
        filename=filename,
        query_string=request.GET.urlencode(),
        genome_build=genome_build,
        vcf_target_system=vcf_target_system,
        included_labs=sorted_labs(export_filter.include_sources or []),
        excluded_labs=sorted_labs(export_filter.exclude_sources or []),
        allele_origin=export_filter.allele_origin,
        since=export_filter.since,
        since_raw=request.GET.get("since") or None,
        rows_per_file=export_settings.rows_per_file,
    )


def view_classification_grouping_export_redirect(request: HttpRequest) -> HttpResponseBase:
    """ Summarises the chosen export settings, then starts the download via a meta refresh to serve_export """
    summary = None
    error = None
    try:
        summary = _export_summary(request)
    except (ValueError, InvalidExportParameter, GenomeBuild.DoesNotExist) as e:
        error = str(e)

    return render(request, "classification/classification_grouping_export_redirect.html", {
        "summary": summary,
        "error": error,
    })


def serve_export(request: HttpRequest) -> HttpResponseBase:
    export_filter = ClassificationGroupingExportFilter.from_request(request)
    export_settings = ClassificationGroupingExportFileSettings.from_request(request)
    classification_export_format = _export_format_from_request(request, export_filter)

    return ClassificationGroupingExportProcess(
        classification_export_format=classification_export_format,
        export_settings=export_settings,
        params=dict(request.GET),
        user=request.user
    ).serve()
