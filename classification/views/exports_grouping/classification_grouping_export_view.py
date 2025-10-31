from django.http import HttpRequest, HttpResponseBase
from django.shortcuts import render
from django.urls import reverse

from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFilter, ClassificationGroupingExportFileSettings
from classification.views.exports_grouping.classification_grouping_export_formatter_csv import \
    ClassificationGroupingExportFormatterCSV
from classification.views.exports_grouping.classification_grouping_export_process import \
    ClassificationGroupingExportProcess
from library.django_utils import get_url_from_view_path
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import GenomeBuild, Lab, Organization


def view_classification_grouping_export(request: HttpRequest) -> HttpResponseBase:
    user_labs = set(lab for lab in Lab.valid_labs_qs(request.user, admin_check=False))
    genome_builds = GenomeBuild.objects.all()
    return render(request, "classification/classification_grouping_export.html", {
        "genome_builds": genome_builds,
        "default_genome_build": GenomeBuildManager.get_current_genome_build(),
        "user_labs": user_labs,
        "all_orgs": list(Organization.objects.filter(active=True).order_by("name")),
        "base_url": get_url_from_view_path(reverse('classification_grouping_export')),
    })


def serve_export(request: HttpRequest) -> HttpResponseBase:
    return ClassificationGroupingExportProcess(
        classification_export_format=ClassificationGroupingExportFormatterCSV(
            ClassificationGroupingExportFilter.from_request(request)
        ),
        export_settings=ClassificationGroupingExportFileSettings()
    ).serve()
