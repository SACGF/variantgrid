from django.db.models.expressions import Subquery
from django.http.request import HttpRequest
from django.shortcuts import render
from rest_framework.response import Response
from termsandconditions.decorators import terms_required

from snpdb.forms import LabSelectForm
from snpdb.models.models_genome import GenomeBuild
from classification.models.classification import Classification, \
    ClassificationModification
from classification.views.classification_datatables import ClassificationDatatableConfig
from classification.views.classification_export_flags import ExportFormatterFlags


@terms_required
def dashboard(request: HttpRequest) -> Response:
    context = {
        "datatable_config": ClassificationDatatableConfig(request),
        "lab_form": LabSelectForm()
    }

    return render(request, 'classification/classification_dashboard.html', context)


def problem_download(request: HttpRequest):
    qs = ClassificationModification.objects.filter(
        is_last_edited=True,
        classification__in=Subquery(
            Classification.filter_for_user(user=request.user).exclude(withdrawn=True).values_list('pk',
                                                                                                  flat=True))
    ).select_related('classification', 'classification__lab')
    exporter = ExportFormatterFlags(
        genome_build=GenomeBuild.grch38(),  # note that genome build for ExportFormatterFlags has no effect
        qs=qs,
        user=request.user
    )
    return exporter.export()
