from crispy_forms.layout import Layout, Field

from django.conf import settings
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.views import View
from django.views.generic import TemplateView

from analysis.models import CandidateSearchRun, CandidateSearchType
from classification.forms import ClassificationAlleleOriginForm, CrossSampleClassificationForm, \
    ClinicalSignificanceForm, ClassificationEvidenceUpdateForm
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


class AbstractNewClassificationCandidateSearchView(View):
    template_name = "classification/candidate_search/abstract_new_classification_candidate_search_view.html"

    def get(self, request):
        return render(request, self.template_name, self.get_context_data())

    def post(self, request):
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
                                                       self._get_candidate_search_type(),
                                                       config_snapshot)
        # Launch a job
        return redirect(csr)

    def _get_layout(self):
        layout_fields = []
        for field_name in ["omim", "hpo", "mondo"]:
            layout_fields.append(Field(field_name, wrapper_class=field_name))
        return Layout(*layout_fields)


    def get_context_data(self):
        # is cached on the request
        user_settings = UserSettingsManager.get_user_settings()

        if settings.CLASSIFICATION_GRID_MULTI_LAB_FILTER:
            lab_form = LabMultiSelectForm()
        else:
            lab_form = LabSelectForm()

        layout = self._get_layout()
        classification_phenotype_form = PhenotypeMultipleSelectForm(prefix="classification")
        classification_phenotype_form.helper.layout = layout

        cs_form = ClinicalSignificanceForm(initial={
            "likely_pathogenic": True,
            "pathogenic": True,
        })

        cs_type = CandidateSearchType(self._get_candidate_search_type())

        return {
            "heading": f"New {cs_type.label} candidate search",
            "button_text": f"Search for {cs_type.label} candidates",
            "gene_form": GeneSymbolForm(prefix="classification"),
            "user_form": UserSelectForm(),
            "lab_form": lab_form,
            "allele_origin_form": ClassificationAlleleOriginForm(),
            "labs": Lab.valid_labs_qs(self.request.user),
            "clinical_significance_form": cs_form,
            "classification_phenotype_form": classification_phenotype_form,
            "user_settings": user_settings,
        }


    def _get_candidate_search_type(self) -> CandidateSearchType:
        raise NotImplementedError()


class NewCrossSampleClassificationCandidateSearchView(AbstractNewClassificationCandidateSearchView):
    template_name = "classification/candidate_search/new_cross_sample_classification_candidate_search.html"

    def get_context_data(self):
        context = super().get_context_data()
        layout = self._get_layout()
        sample_phenotype_form = PhenotypeMultipleSelectForm(prefix="sample")
        sample_phenotype_form.helper.layout = layout
        context["sample_phenotype_form"] = sample_phenotype_form
        context["sample_classification_form"] = CrossSampleClassificationForm()

        samples_qs = Sample.filter_for_user(self.request.user)
        context["num_visible_samples"] = samples_qs.count()
        num_sample_gene_lists = SampleGeneList.objects.filter(sample__in=samples_qs).count()
        context["num_sample_gene_lists"] = num_sample_gene_lists
        if num_sample_gene_lists:
            sample_gene_form = GeneSymbolForm(prefix="sample")
            sample_gene_form.fields['gene_symbol'].label = "Sample Gene List Gene Symbol"
            context["sample_gene_form"] = sample_gene_form

        return context

    def _get_candidate_search_type(self) -> CandidateSearchType:
        return CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION


class NewClassificationEvidenceUpdateCandidateSearchView(AbstractNewClassificationCandidateSearchView):
    template_name = "classification/candidate_search/new_classification_evidence_update_candidate_search.html"

    def get_context_data(self):
        context = super().get_context_data()
        context["classification_evidence_update_form"] = ClassificationEvidenceUpdateForm()
        return context

    def _get_candidate_search_type(self) -> CandidateSearchType:
        return CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE

