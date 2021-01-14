import urllib

from django.contrib import messages
from django.views.generic import TemplateView
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake
from ontology.ontology_matching import OntologyMatching


class SearchMondoText(APIView):

    def get(self, request, **kwargs) -> Response:

        search_term = request.GET.get('search_term') or ''
        gene_symbol = request.GET.get('gene_symbol')

        urllib.parse.quote(search_term).replace('/', '%252F')  # a regular escape / gets confused for a URL divider
        selected = [term.strip() for term in (request.GET.get('selected') or '').split(",") if term.strip()]

        ontology_matches = OntologyMatching.from_search(search_text=search_term, gene_symbol=gene_symbol, selected=selected, server_search=True)

        return Response(status=HTTP_200_OK, data=ontology_matches.as_json())


class OntologyTermView(TemplateView):

    template_name = "ontology/ontology_term.html"

    def get_context_data(self, **kwargs):
        term_id = self.kwargs.get("term")
        term = OntologyTerm.get_or_stub(term_id)
        if not term.is_stub:
            gene_relationships = None
            if term.ontology_service != OntologyService.HGNC:
                gene_relationships = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC)
            return {
                "term": term,
                "gene_relationships": gene_relationships,
                "relationships": OntologyTermRelation.relations_of(term)
            }
        else:
            messages.add_message(self.request, messages.ERROR, "This term is not stored in our database")
            return {"term": term}
