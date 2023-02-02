import urllib

from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from genes.models import GeneListGeneSymbol, create_fake_gene_list
from genes.serializers import GeneListGeneSymbolSerializer
from library.constants import WEEK_SECS
from ontology.models import OntologyTerm, GeneDiseaseClassification, OntologyVersion
from ontology.ontology_matching import OntologyMatching
from ontology.serializers import OntologyTermRelationSerializer


class SearchMondoText(APIView):
    def get(self, request, **kwargs) -> Response:

        search_term = request.GET.get('search_term') or ''
        gene_symbol = request.GET.get('gene_symbol')

        urllib.parse.quote(search_term).replace('/', '%252F')  # a regular escape / gets confused for a URL divider
        selected = [term.strip() for term in (request.GET.get('selected') or '').split(",") if term.strip()]

        ontology_matches = OntologyMatching.from_search(search_text=search_term, gene_symbol=gene_symbol, selected=selected)

        return Response(status=HTTP_200_OK, data=ontology_matches.as_json())


class OntologyTermGeneListView(APIView):
    def get(self, request, *args, **kwargs):
        term_slug = self.kwargs['term']
        ontology_version = OntologyVersion.latest()
        ontology_term = OntologyTerm.get_from_slug(term_slug)

        name = str(ontology_term)
        gene_list = create_fake_gene_list(name=name, user=None)
        gene_list_genes = []
        for gene_symbol in ontology_version.gene_symbols_for_terms([ontology_term.pk]):
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


@method_decorator(cache_page(WEEK_SECS), name='get')
class GeneDiseaseRelationshipView(APIView):
    def get(self, request, *args, **kwargs):
        data = []
        ontology_version = OntologyVersion.latest()
        for otr in ontology_version.gene_disease_relations(self.kwargs['gene_symbol'],
                                                           min_classification=GeneDiseaseClassification.DISPUTED):
            data.append(OntologyTermRelationSerializer(otr).data)
        return Response(data)
