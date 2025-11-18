from functools import cached_property

from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404

from analysis.models import CandidateSearchRun, Candidate
from classification.views.views import CreateClassificationForVariantView
from snpdb.forms import SampleChoiceForm
from snpdb.models import GenomeBuild, Sample


def view_candidate_search_run(request, pk) -> HttpResponse:
    """ This is the generic results page """
    candidate_search_run = get_object_or_404(CandidateSearchRun, pk=pk)
    # Permission check??
    context = {
        "candidate_search_run": candidate_search_run,
    }
    return render(request, 'analysis/candidate_search/view_candidate_search_run.html', context)


def reanalyis(request):
    # Another name for this would be analysis candidate search
    context = {}
    return render(request, 'analysis/candidate_search/reanalysis.html', context)


def new_reanalyis_candidate_search(request):
    context = {
    }
    return render(request, 'analysis/candidate_search/new_reanalysis_candidate_search.html', context)




class CreateClassificationForCandidateView(CreateClassificationForVariantView):
    template_name = "analysis/candidate_search/create_classification_for_candidate.html"

    def _get_variant(self):
        return self.candidate.variant

    def _get_genome_build(self) -> GenomeBuild:
        return self.candidate.sample.genome_build

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
        candidate_id = self.kwargs["candidate_id"]
        c = Candidate.objects.get(pk=candidate_id)
        c.search_run.check_can_view(self.request.user)  # Permission check
        return c

    def get_context_data(self, *args, **kwargs):
        context = super().get_context_data(*args, **kwargs)
        context["candidate"] = self.candidate
        return context
