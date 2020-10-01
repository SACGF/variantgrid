from django.db.models.expressions import Subquery
from django.http.request import HttpRequest
from django.shortcuts import render
from rest_framework.response import Response
from termsandconditions.decorators import terms_required

from snpdb.forms import LabSelectForm
from snpdb.models.models_genome import GenomeBuild
from classification.models.variant_classification import VariantClassification, \
    VariantClassificationModification
from classification.views.variant_classification_datatables import VariantClassificationDatatableConfig
from classification.views.variant_classification_export_flags import ExportFormatterFlags


@terms_required
def dashboard(request: HttpRequest) -> Response:
    context = {
        "datatable_config": VariantClassificationDatatableConfig(request),
        "lab_form": LabSelectForm()
    }

    return render(request, 'classification/variant_classification_dashboard.html', context)


def problem_download(request: HttpRequest):
    qs = VariantClassificationModification.objects.filter(
        is_last_edited=True,
        variant_classification__in=Subquery(
            VariantClassification.filter_for_user(user=request.user).exclude(withdrawn=True).values_list('pk',
                                                                                                         flat=True))
    ).select_related('variant_classification', 'variant_classification__lab')
    exporter = ExportFormatterFlags(
        genome_build=GenomeBuild.grch38(),  # note that genome build for ExportFormatterFlags has no effect
        qs=qs,
        user=request.user
    )
    return exporter.export()
