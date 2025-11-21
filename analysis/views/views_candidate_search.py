import logging
from functools import cached_property

from django.core.exceptions import PermissionDenied
from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from django.views.decorators.http import require_POST
from django.views.generic.base import TemplateView

from analysis.forms import ReanalysisCandidateSearchForm, AnalysisFilterForm
from analysis.models import CandidateSearchRun, Candidate, CandidateStatus, CandidateSearchType
from classification.views.views import CreateClassificationForVariantView, create_classification_object
from snpdb.forms import SampleChoiceForm
from snpdb.models import GenomeBuild


def view_candidate_search_run(request, pk) -> HttpResponse:
    """ This is the generic results page """
    candidate_search_run = get_object_or_404(CandidateSearchRun, pk=pk)
    # Permission check??
    context = {
        "candidate_search_run": candidate_search_run,
    }
    return render(request, 'analysis/candidate_search/view_candidate_search_run.html', context)


class AbstractCandidateSearchView(TemplateView):
    def _get_search_types(self) -> list[CandidateSearchType]:
        raise NotImplementedError()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["search_types"] = self._get_search_types()
        return context


class ReanalyisCandidateSearchView(AbstractCandidateSearchView):
    template_name = "analysis/candidate_search/reanalysis_candidate_search.html"
    def _get_search_types(self) -> list[CandidateSearchType]:
        return [CandidateSearchType.REANALYSIS_NEW_ANNOTATION]


def new_reanalyis_candidate_search(request):
    search_type = CandidateSearchType.REANALYSIS_NEW_ANNOTATION
    context = {
        "heading": f"Search for new {search_type.label}",
        "methods": search_type.get_methods(),
        "button_text": "Search for reanalysis candidates",
        "analyses_filter_form": AnalysisFilterForm(),
        "reanalysis_candidate_search_form": ReanalysisCandidateSearchForm(),
    }
    return render(request, 'analysis/candidate_search/new_reanalysis_candidate_search.html', context)


class CreateClassificationForCandidateView(CreateClassificationForVariantView):
    """ This is a wrapper around Classification """
    # I am not sure whether this belongs in Anaylsis or just Classification - will see if we ever want to classify
    # something in analysis candidate results
    template_name = "analysis/candidate_search/create_classification_for_candidate.html"

    def _get_variant(self):
        return self.candidate.variant

    def _get_genome_build(self) -> GenomeBuild:
        return self.candidate.sample.genome_build

    def _get_form_post_url(self) -> str:
        return reverse("create_classification_for_candidate", kwargs={"candidate_id": self.candidate.pk})

    def _get_sample_form(self):
        # Use candidate sample - otherwise fall back on default (all samples in DB visible to user)
        if self.candidate.sample:
            form = SampleChoiceForm(initial={"sample": self.candidate.sample})
            form.fields['sample'].disabled = True
        else:
            form = super()._get_sample_form()
        return form

    @cached_property
    def candidate(self):
        return Candidate.get_permission_check(pk=self.kwargs["candidate_id"], user=self.request.user)

    def get_context_data(self, *args, **kwargs):
        context = super().get_context_data(*args, **kwargs)
        context["candidate"] = self.candidate
        return context


@require_POST
def create_classification_for_candidate(request, candidate_id):
    """ Performs classification then closes candidate (if have write permissions) """
    classification = create_classification_object(request)
    try:
        candidate = Candidate.get_permission_check(pk=candidate_id, user=request.user, write=True)
        candidate.reviewer = request.user
        candidate.notes = f"Created classification: {classification.pk}"
        candidate.status = CandidateStatus.RESOLVED
        candidate.save()
    except PermissionDenied:
        # This doesn't actually matter much - just can't resolve it
        logging.warning("%s does not have permission to write to candidate=%s", request.user, candidate_id)
    return redirect(classification.get_edit_url())
