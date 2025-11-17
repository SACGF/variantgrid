from crispy_forms.layout import Layout, Field

from django.conf import settings
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect

from analysis.models import CandidateSearchRun, CandidateSearchType
from classification.forms import ClassificationAlleleOriginForm, SampleClassificationForm, ClinicalSignificanceForm
from genes.forms import GeneSymbolForm
from genes.models import SampleGeneList
from ontology.forms import PhenotypeMultipleSelectForm
from snpdb.forms import UserSelectForm, LabSelectForm, LabMultiSelectForm
from snpdb.models import Lab, Sample
from snpdb.user_settings_manager import UserSettingsManager


def view_classification_candidate_search(request, pk) -> HttpResponse:
    classification_candidate_search_run = get_object_or_404(CandidateSearchRun, pk=pk)
    # Permission check??
    context = {
        "classification_candidate_search_run": classification_candidate_search_run,
    }
    return render(request, 'classification/candidate_search/classification_candidate_search.html', context)



def classification_candidate_search(request) -> HttpResponse:
    context = {}
    return render(request, 'classification/candidate_search/classification_candidate_search.html', context)


def new_cross_sample_classification_candidate_search(request) -> HttpResponse:
    """ TODO: This could be a class then use inheritance to handle grid being in common?? """

    if request.method == "POST":
        post_ignore_contains = [
            "csrfmiddlewaretoken",
            "datatable_length",
        ]
        config_snapshot = {}
        for k, v in request.POST.items():
            ignore = False
            for ignore_str in post_ignore_contains:
                if ignore := ignore_str in k:
                    break
            if not ignore:
                config_snapshot[k] = v

        csr = CandidateSearchRun.create_and_launch_job(request.user,
                                                       CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION,
                                                       config_snapshot)
        # Launch a job
        return redirect(csr)
    else:
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

        cs_form = ClinicalSignificanceForm(initial={
            "likely_pathogenic": True,
            "pathogenic": True,
        })

        context = {
            "sample_phenotype_form": sample_phenotype_form,
            "gene_form": GeneSymbolForm(prefix="classification"),
            "user_form": UserSelectForm(),
            "lab_form": lab_form,
            "allele_origin_form": ClassificationAlleleOriginForm(),
            "labs": Lab.valid_labs_qs(request.user),
            "clinical_significance_form": cs_form,
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
        return render(request, 'classification/candidate_search/new_cross_sample_classification_candidate_search.html', context)


def new_classification_evidence_update_candidate_search(request) -> HttpResponse:
    context = {}
    return render(request, 'classification/candidate_search/new_classification_evidence_update_candidate_search.html', context)
