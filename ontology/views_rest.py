import urllib

from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from genes.models import FakeGeneList, GeneListGeneSymbol
from genes.serializers import GeneListGeneSymbolSerializer
from ontology.models import OntologyTerm, OntologySnake
from ontology.ontology_matching import OntologyMatching


class SearchMondoText(APIView):
    def get(self, request, **kwargs) -> Response:

        search_term = request.GET.get('search_term') or ''
        gene_symbol = request.GET.get('gene_symbol')

        urllib.parse.quote(search_term).replace('/', '%252F')  # a regular escape / gets confused for a URL divider
        selected = [term.strip() for term in (request.GET.get('selected') or '').split(",") if term.strip()]

        ontology_matches = OntologyMatching.from_search(search_text=search_term, gene_symbol=gene_symbol, selected=selected, server_search=True)

        return Response(status=HTTP_200_OK, data=ontology_matches.as_json())


class OntologyTermGeneListView(APIView):
    def get(self, request, *args, **kwargs):
        term_slug = self.kwargs['term']
        ontology_term = OntologyTerm.get_from_slug(term_slug)

        name = str(ontology_term)
        gene_list = FakeGeneList(name=name, user=None)
        gene_list_genes = []
        for gene_symbol in OntologySnake.gene_symbols_for_term(ontology_term):
            glg = GeneListGeneSymbol(gene_list=gene_list, original_name=gene_symbol, gene_symbol=gene_symbol)
            gene_list_genes.append(glg)

        genelistgenesymbol_set = GeneListGeneSymbolSerializer(gene_list_genes, many=True).data
        ontology_service_name = ontology_term.get_ontology_service_display()

        data = {"pk": f"ontology-{term_slug}",
                "category": {"name": ontology_service_name,
                             "icon_css_class": f"{ontology_service_name}-icon"},
                "name": name,
                "import_status": "S",
                "genelistgenesymbol_set": genelistgenesymbol_set,
                "can_write": False}
        return Response(data)
