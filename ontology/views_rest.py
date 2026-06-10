import re

from django.http import Http404
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from genes.models import GeneListGeneSymbol, create_fake_gene_list
from genes.serializers import GeneListGeneSymbolSerializer
from library.constants import WEEK_SECS
from ontology.models import (
    ONTOLOGY_RELATIONSHIP_MEDIUM_QUALITY_FILTER,
    OntologyTerm,
    OntologyVersion,
)
from ontology.ontology_matching import OntologyMatching
from ontology.serializers import OntologyTermRelationSerializer

_GENE_SYMBOL_RE = re.compile(r'^[A-Za-z0-9\-]+$')


class SearchMondoText(APIView):
    """ Searches MONDO ontology terms by free text, optionally scoped to a gene symbol """

    @extend_schema(
        summary="Search MONDO ontology terms by text, optionally filtered by gene symbol",
        parameters=[
            OpenApiParameter("search_term", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Free text to search MONDO terms for"),
            OpenApiParameter("gene_symbol", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Gene symbol to scope the search to"),
            OpenApiParameter("selected", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Comma-separated list of already-selected ontology term IDs"),
        ],
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, **kwargs) -> Response:

        search_term = request.GET.get('search_term') or ''
        gene_symbol = request.GET.get('gene_symbol')

        selected = [term.strip() for term in (request.GET.get('selected') or '').split(",") if term.strip()]

        ontology_matches = OntologyMatching.from_search(search_text=search_term, gene_symbol=gene_symbol, selected=selected)

        return Response(status=HTTP_200_OK, data=ontology_matches.as_json())


class OntologyTermGeneListView(APIView):
    """ Returns a fake gene list built from the gene symbols associated with an ontology term """

    @extend_schema(
        summary="Retrieve a gene list of gene symbols associated with an ontology term",
        responses=OpenApiTypes.OBJECT,
    )
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
    """ Returns gene/disease relationships (medium quality or better) for a gene symbol """

    @extend_schema(
        summary="List gene/disease ontology relationships for a gene symbol",
        parameters=[
            OpenApiParameter("gene_symbol", OpenApiTypes.STR, OpenApiParameter.PATH,
                             description="Gene symbol, e.g. BRCA1"),
        ],
        responses=OntologyTermRelationSerializer(many=True),
    )
    def get(self, request, *args, **kwargs):
        gene_symbol = self.kwargs['gene_symbol']
        if not _GENE_SYMBOL_RE.match(gene_symbol):
            raise Http404
        data = []
        ontology_version = OntologyVersion.latest()
        for otr in ontology_version.gene_disease_relations(gene_symbol, quality_filter=ONTOLOGY_RELATIONSHIP_MEDIUM_QUALITY_FILTER):
            data.append(OntologyTermRelationSerializer(otr).data)
        return Response(data)
