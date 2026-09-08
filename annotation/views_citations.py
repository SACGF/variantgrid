from django.http.response import JsonResponse
from django.shortcuts import render
from django.template.loader import render_to_string

from annotation.models import Citation
from annotation.models.models_citations import CitationFetchRequest
from library.utils import first


def view_citation(request, citation_id: str):
    """
    A dedicated page for a single citation, handy for testing.
    :param request: The Request
    :param citation_id: The full ID of the citation, e.g. PMID:4353345
    :return: A stand alone page with the loaded citation details
    """
    # important, we don't actually retrieve the citation here
    return render(request, "annotation/citation.html", {"citation": first(CitationFetchRequest.get_unfetched_citations([citation_id]))})


def view_citation_detail(request, citation_id: str):
    """
    A dedicated page for a single citation, handy for testing
    :param request: The Request
    :param citation_id: The full ID of the citation, e.g. PMID:4353345
    :return: A stand alone page with the loaded citation details
    """
    return render(request, "annotation/citation_detail.html", {"citation": CitationFetchRequest.fetch_all_now([citation_id]).first_citation})


def simple_citation_html(cd: Citation) -> str:
    return render_to_string('annotation/citation_simple.html', {"citation": cd}).replace('\n', '').strip()


def citations_json(request, citations_ids_list):
    """
    Request JSON of citations, accepts either variant grid internal citation ids, or PubMed:123456
    """
    citation_ids = citations_ids_list.split("/")
    return JsonResponse({"citations": CitationFetchRequest.fetch_all_now(citation_ids).to_json()})
