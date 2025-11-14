from django.http.response import HttpResponse
from django.shortcuts import render, get_object_or_404

from analysis.models import CandidateSearchRun


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


