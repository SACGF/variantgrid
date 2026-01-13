from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Field
from django.conf import settings
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404

from analysis.forms import SampleCandidatesSearchForm
from analysis.models import CandidateSearchRun, CandidateSearchType
from analysis.views.views_candidate_search import AbstractCandidateSearchView, AbstractNewCandidateSearchView
from classification.forms import ClassificationAlleleOriginForm, ClinicalSignificanceForm, \
    ClassificationEvidenceUpdateForm
from genes.forms import GeneSymbolForm
from genes.models import SampleGeneList
from ontology.forms import PhenotypeMultipleSelectForm
from snpdb.forms import UserSelectForm, LabSelectForm, LabMultiSelectForm, ProjectChoiceForm, VCFChoiceForm, \
    SampleMultiForm
from snpdb.models import Lab, Sample
from snpdb.user_settings_manager import UserSettingsManager


def view_classification_candidate_search(request, pk) -> HttpResponse:
    classification_candidate_search_run = get_object_or_404(CandidateSearchRun, pk=pk)
    # Permission check??
    context = {
        "classification_candidate_search_run": classification_candidate_search_run,
    }
    return render(request, 'classification/candidate_search/classification_candidate_search.html', context)


class ReanalyisCandidateSearchView(AbstractCandidateSearchView):
    template_name = 'classification/candidate_search/classification_candidate_search.html'
    def _get_search_types(self) -> list[CandidateSearchType]:
        return [CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION, CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE]


class AbstractNewClassificationCandidateSearchView(AbstractNewCandidateSearchView):
    template_name = "classification/candidate_search/abstract_new_classification_candidate_search_view.html"

    def _get_layout(self):
        layout_fields = []
        for field_name in ["omim", "hpo", "mondo"]:
            layout_fields.append(Field(field_name, wrapper_class=field_name))
        return Layout(*layout_fields)

    def get_context_data(self):
        context = super().get_context_data()
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

        context.update({
            "gene_form": GeneSymbolForm(prefix="classification"),
            "user_form": UserSelectForm(),
            "lab_form": lab_form,
            "allele_origin_form": ClassificationAlleleOriginForm(),
            "labs": Lab.valid_labs_qs(self.request.user),
            "clinical_significance_form": cs_form,
            "classification_phenotype_form": classification_phenotype_form,
            "user_settings": user_settings,
        })
        return context


class NewCrossSampleClassificationCandidateSearchView(AbstractNewClassificationCandidateSearchView):
    template_name = "classification/candidate_search/new_cross_sample_classification_candidate_search.html"

    def get_context_data(self):
        context = super().get_context_data()
        layout = self._get_layout()

        no_form_helper = FormHelper()
        no_form_helper.form_tag = False
        sample_candidate_search_form = SampleCandidatesSearchForm()
        sample_candidate_search_form.helper = no_form_helper
        project_choice_form = ProjectChoiceForm()
        project_choice_form.helper = no_form_helper
        vcf_choice_form = VCFChoiceForm()
        vcf_choice_form.helper = no_form_helper
        sample_multiple_choice_form = SampleMultiForm()
        sample_multiple_choice_form.helper = no_form_helper

        sample_phenotype_form = PhenotypeMultipleSelectForm(prefix="sample")
        sample_phenotype_form.helper.layout = layout
        context["sample_phenotype_form"] = sample_phenotype_form
        context["sample_candidate_search_form"] = sample_candidate_search_form
        context["project_choice_form"] = project_choice_form
        context["vcf_choice_form"] = vcf_choice_form
        context["sample_multiple_choice_form"] = sample_multiple_choice_form

        samples_qs = Sample.filter_for_user(self.request.user)
        context["num_visible_samples"] = samples_qs.count()
        num_sample_gene_lists = SampleGeneList.objects.filter(sample__in=samples_qs).count()
        context["num_sample_gene_lists"] = num_sample_gene_lists
        if num_sample_gene_lists:
            sample_gene_form = GeneSymbolForm(prefix="sample")
            sample_gene_form.helper = no_form_helper
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
