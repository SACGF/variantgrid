from django.db.models.expressions import Subquery
from django.http.request import HttpRequest
from django.shortcuts import render
from rest_framework.response import Response
from termsandconditions.decorators import terms_required

from classification.models import classification_flag_types
from flags.models import FlagCollection
from snpdb.forms import LabSelectForm
from snpdb.models import Lab
from snpdb.models.models_genome import GenomeBuild
from classification.models.classification import Classification, \
    ClassificationModification
from classification.views.classification_datatables import ClassificationDatatableConfig
from classification.views.classification_export_flags import ExportFormatterFlags


@terms_required
def dashboard(request: HttpRequest) -> Response:
    classification_counts = dict()

    labs = list(Lab.valid_labs_qs(request.user, admin_check=True))
    labs.sort()

    vcqs_user = Classification.filter_for_user(request.user).filter(lab__in=labs)
    vcqs = vcqs_user.filter(withdrawn=False)

    discordant = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
        classification_flag_types.discordant
    ])
    variant_matching = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
        classification_flag_types.matching_variant_flag,
        classification_flag_types.matching_variant_warning_flag,
        classification_flag_types.transcript_version_change_flag
    ])
    comments = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
        classification_flag_types.suggestion,
        classification_flag_types.internal_review,
        classification_flag_types.significance_change
    ])
    unshared = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[
        classification_flag_types.unshared_flag
    ])
    withdrawn = FlagCollection.filter_for_open_flags(qs=vcqs_user.filter(withdrawn=True), flag_types=[
        classification_flag_types.classification_withdrawn
    ])

    issue_counts = {
        "classifications": {
            "discordant": discordant.count(),
            "variant_matching": variant_matching.count(),
            "comments": comments.count(),
            "unshared": unshared.count(),
            "withdrawn": withdrawn.count()
        }
    }

    context = {
        "datatable_config": ClassificationDatatableConfig(request),
        "labs": labs,
        "counts": issue_counts
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
