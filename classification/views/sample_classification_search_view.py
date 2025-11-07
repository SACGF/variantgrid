from crispy_forms.layout import Layout, Field

from django.conf import settings
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.shortcuts import render

from classification.forms import ClassificationAlleleOriginForm, SampleClassificationForm
from genes.forms import GeneSymbolForm
from genes.models import SampleGeneList
from ontology.forms import PhenotypeMultipleSelectForm
from snpdb.forms import UserSelectForm, LabSelectForm, LabMultiSelectForm
from snpdb.models import Lab, Sample
from snpdb.user_settings_manager import UserSettingsManager


def sample_classification_search(request) -> HttpResponse:
    # A lot of this code is copy/pasted from classifications - classification listing page
    # TODO: We should go and refactor this into a tag later

    # is cached on the request
    user_settings = UserSettingsManager.get_user_settings()

    if settings.CLASSIFICATION_GRID_MULTI_LAB_FILTER:
        lab_form = LabMultiSelectForm()
    else:
        lab_form = LabSelectForm()

    sample_gene_form = GeneSymbolForm(prefix="sample")
    sample_gene_form.fields['gene_symbol'].label = "Sample Gene List Gene Symbol"

    layout_fields = []
    for field_name in ["omim", "hpo", "mondo"]:
        layout_fields.append(Field(field_name, wrapper_class=field_name))
    layout = Layout(*layout_fields)

    sample_phenotype_form = PhenotypeMultipleSelectForm(prefix="sample")
    sample_phenotype_form.helper.layout = layout

    classification_phenotype_form = PhenotypeMultipleSelectForm(prefix="classification")
    classification_phenotype_form.helper.layout = layout

    context = {
        "sample_phenotype_form": sample_phenotype_form,
        "gene_form": GeneSymbolForm(prefix="classification"),
        "user_form": UserSelectForm(),
        "lab_form": lab_form,
        "allele_origin_form": ClassificationAlleleOriginForm(),
        "labs": Lab.valid_labs_qs(request.user),
        "classification_phenotype_form": classification_phenotype_form,
        "sample_classification_form": SampleClassificationForm(),
        "user_settings": user_settings,
    }

    samples_qs = Sample.filter_for_user(request.user)
    context["num_visible_samples"] = samples_qs.count()
    num_sample_gene_lists = SampleGeneList.objects.filter(sample__in=samples_qs).count()
    context["num_sample_gene_lists"] = num_sample_gene_lists
    if num_sample_gene_lists:
        context["sample_gene_form"] = sample_gene_form
    return render(request, 'classification/sample_classification_search.html', context)


def sample_classification_search_results(request: HttpRequest) -> HttpResponse:
    # Sample filters
    request.GET.get("sample_gene_symbol")
    request.GET.get("sample_ontology_term_id")

    # Classification filters
    request.GET.get("classification_allele_origin")
    request.GET.get("classification_gene_symbol")
    request.GET.get("classification_id_filter")
    request.GET.get("classification_lab")
    request.GET.get("classification_ontology_term_id")
    request.GET.get("classification_user")
    # Search
    request.GET.get("search_max_results")
    request.GET.get("search_max_samples")
    context = {

    }
    return render(request, 'classification/sample_classification_search_results.html', context)
