from django.conf import settings
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.shortcuts import render

from classification.forms import ClassificationAlleleOriginForm
from genes.forms import GeneSymbolForm
from snpdb.forms import UserSelectForm, LabSelectForm, LabMultiSelectForm
from snpdb.models import Lab
from snpdb.user_settings_manager import UserSettingsManager


def sample_classification_search(request: HttpRequest) -> HttpResponse:
    # A lot of this code is copy/pasted from classifications - classification listing page
    # TODO: We should go and refactor this into a tag later

    # is cached on the request
    user_settings = UserSettingsManager.get_user_settings()

    if settings.CLASSIFICATION_GRID_MULTI_LAB_FILTER:
        lab_form = LabMultiSelectForm()
    else:
        lab_form = LabSelectForm()

    context = {
        "gene_form": GeneSymbolForm(),
        "user_form": UserSelectForm(),
        "lab_form": lab_form,
        "allele_origin_form": ClassificationAlleleOriginForm(),
        "labs": Lab.valid_labs_qs(request.user),
        "genome_build": user_settings.default_genome_build,
        "user_settings": user_settings,
    }
    return render(request, 'classification/sample_classification_search.html', context)

