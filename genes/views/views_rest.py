import json
import logging
from collections import defaultdict

from django.http.response import Http404
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema, OpenApiParameter
from rest_framework import permissions
from rest_framework.generics import ListAPIView
from rest_framework.generics import RetrieveAPIView
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from genes.gene_matching import GeneSymbolMatcher
from genes.models import Gene, GeneInfo, GeneList, GeneAnnotationRelease, GeneSymbol, \
    ReleaseGeneSymbolGene, PanelAppServer, SampleGeneList, ActiveSampleGeneList, create_fake_gene_list
from genes.panel_app import get_panel_app_panel_as_gene_list_json
from genes.panel_app import get_panel_app_results_by_gene_symbol_json
from genes.serializers import GeneDetailSerializer, GeneInfoSerializer, GeneListGeneSymbolSerializer, \
    GeneListSerializer, GeneAnnotationReleaseSerializer, GeneSymbolDetailSerializer, SampleGeneListSerializer
from library.constants import HOUR_SECS, WEEK_SECS
from library.django_utils.django_rest_utils import MultipleFieldLookupMixin
from library.guardian_utils import DjangoPermission
from snpdb.models import GenomeBuild
from snpdb.models.models_enums import ImportStatus

log = logging.getLogger(__name__)


def is_owner_or_has_permission_factory(django_permission):

    class IsOwnerOrHasReadPermission(permissions.BasePermission):

        def has_object_permission(self, request, view, obj):
            is_owner = obj.user == request.user
            permission = DjangoPermission.perm(obj, django_permission)
            has_permission = request.user.has_perm(permission, obj)
            return is_owner or has_permission or request.user.is_superuser

    return IsOwnerOrHasReadPermission


WriteGeneListPermission = is_owner_or_has_permission_factory(DjangoPermission.WRITE)


class PanelAppGeneListView(APIView):
    """ Tunnels through to panel app (can't make cross site requests) """

    @extend_schema(
        summary="Retrieve a PanelApp panel converted to gene list JSON",
        parameters=[OpenApiParameter("id", OpenApiTypes.INT, OpenApiParameter.PATH,
                                     description="PanelApp panel ID")],
        responses=OpenApiTypes.OBJECT,
    )
    @method_decorator(cache_page(HOUR_SECS))
    def get(self, request, *args, **kwargs):
        panel_app_id = self.kwargs['pk']
        data = get_panel_app_panel_as_gene_list_json(panel_app_id)
        return Response(data)


class GeneListView(APIView):
    """ Retrieve a gene list (with gene symbols) the user has permission to view """

    @extend_schema(
        summary="Retrieve a gene list by ID",
        parameters=[OpenApiParameter("id", OpenApiTypes.INT, OpenApiParameter.PATH,
                                     description="GeneList ID")],
        responses=GeneListSerializer,
    )
    def get(self, request, *args, **kwargs):
        gene_list_id = self.kwargs['pk']
        gl = GeneList.get_for_user(request.user, gene_list_id)
        serializer = GeneListSerializer(gl, context={"request": request})
        data = serializer.data

        return Response(data)


class NamedGeneListView(MultipleFieldLookupMixin, RetrieveAPIView):
    serializer_class = GeneListSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    lookup_fields = ('category__name', 'name')

    def get_queryset(self):
        #category = self.kwargs['category']
        #name = self.kwargs['name']
        return GeneList.filter_for_user(self.request.user)  # .filter(category__name=category, name=name)


class ModifyGeneListView(APIView):
    """ Add/remove gene symbols from an existing gene list the user can write to """
    permission_classes = (permissions.IsAuthenticated, WriteGeneListPermission)

    def get_gene_list_modifications(self, data):
        modifications = json.loads(data["modifications"])
        gene_additions = modifications.get('add')
        gene_deletions = modifications.get('delete')
        return gene_additions, gene_deletions

    @extend_schema(
        summary="Add and/or remove gene symbols from a gene list",
        request=OpenApiTypes.OBJECT,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, **kwargs):
        gene_list_id = kwargs.pop("pk")
        gene_list = GeneList.get_for_user(request.user, gene_list_id)
        self.check_object_permissions(request, gene_list)
        gene_additions, gene_deletions = self.get_gene_list_modifications(request.data)

        modification_info = f"Added manually by {request.user} on {timezone.now()}"
        gene_additions_modification_info = {gene: modification_info for gene in gene_additions}
        num_added, num_deleted = gene_list.add_and_remove_gene_symbols(gene_additions, gene_deletions,
                                                                       gene_additions_modification_info)
        return Response(status=HTTP_200_OK, data={"num_added": num_added, "num_deleted": num_deleted})


class PanelAppGeneEvidenceView(APIView):
    """ Need to tunnel through due to Cross site requests """

    @extend_schema(
        summary="Retrieve PanelApp evidence for a gene symbol from a PanelApp server",
        parameters=[OpenApiParameter("gene_symbol", OpenApiTypes.STR, OpenApiParameter.PATH,
                                     description="Gene symbol to look up")],
        responses=OpenApiTypes.OBJECT,
    )
    @method_decorator(cache_page(WEEK_SECS))
    def get(self, request, *args, **kwargs):
        server_id = self.kwargs['server_id']
        gene_symbol = self.kwargs['gene_symbol']
        server = PanelAppServer.objects.get(pk=server_id)
        data = get_panel_app_results_by_gene_symbol_json(server, gene_symbol)
        if data is None:
            raise Http404(f"PanelApp has no results for '{gene_symbol}'")
        return Response(data)


class CreateGeneListView(APIView):
    """ Create a new gene list from a name and list of gene symbols """

    @extend_schema(
        summary="Create a gene list from gene symbols",
        request=OpenApiTypes.OBJECT,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request):
        data = request.data
        name = data["name"]
        gene_symbols = data["gene_symbols"]
        modification_info = data.get("modification_info", "Created via API")

        gene_list = GeneList.objects.create(name=name, user=request.user)
        try:
            gene_matcher = GeneSymbolMatcher()
            gene_matcher.create_gene_list_gene_symbols(gene_list, gene_symbols, modification_info)
            import_status = ImportStatus.SUCCESS
        except Exception:
            log.exception("Error creating gene list %d for user %s", gene_list.pk, request.user)
            gene_list.error_message = "An error occurred while importing the gene list."
            import_status = ImportStatus.ERROR

        gene_list.import_status = import_status
        gene_list.save()

        return Response({"pk": gene_list.pk})


class TextToGeneListView(APIView):
    """ Text to gene list (doesn't actually save). Used to e.g. check gene names """

    @extend_schema(
        summary="Match free text gene symbols as an unsaved gene list (name checking)",
        request=OpenApiTypes.OBJECT,
        responses=GeneListSerializer,
    )
    def post(self, request, *args, **kwargs):
        # Needed to be post as ~500 genes exceeded GET limit of ~4k
        name = self.request.data['name']
        gene_list_text = self.request.data['gene_list_text']

        gene_list = create_fake_gene_list(name=name, user=request.user)
        context = {
            "request": request,
            "exclude_fields": ["genelistgenesymbol_set"],
        }
        serializer = GeneListSerializer(gene_list,
                                        context=context)
        data = serializer.data

        if gene_list_text:
            gene_matcher = GeneSymbolMatcher()
            gene_list_gene_symbols = gene_matcher.create_gene_list_gene_symbols_from_text(gene_list, gene_list_text,
                                                                                          save=False)
            data["genelistgenesymbol_set"] = GeneListGeneSymbolSerializer(gene_list_gene_symbols, many=True).data
        return Response(data)


class GeneAnnotationReleaseView(RetrieveAPIView):
    serializer_class = GeneAnnotationReleaseSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    lookup_fields = ('pk',)

    def get_queryset(self):
        return GeneAnnotationRelease.objects.all()


@extend_schema(
    parameters=[OpenApiParameter("gene_symbol", OpenApiTypes.STR, OpenApiParameter.PATH,
                                 description="Gene symbol to retrieve GeneInfo for")],
)
class GeneInfoView(ListAPIView):
    """ List GeneInfo records (e.g. tags/icons) for a gene symbol """
    serializer_class = GeneInfoSerializer

    def get_queryset(self):
        gene_symbol = self.kwargs['gene_symbol']
        return GeneInfo.get_for_gene_symbol(gene_symbol)


class BatchGeneInfoView(APIView):
    """ Needs to be a post as we can send a large number of genes
        returns {gene_symbol : [gene_info_dict1, gene_info_dict2]} """

    @extend_schema(
        summary="Retrieve GeneInfo records for a batch of gene symbols",
        request=OpenApiTypes.OBJECT,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request):
        gene_symbols_json = request.data["gene_symbols_json"]
        gene_symbols = json.loads(gene_symbols_json)
        gene_symbol_path = "gene_list__genelistgenesymbol__gene_symbol"

        qs = GeneInfo.objects.filter(**{gene_symbol_path + "__in": gene_symbols}).distinct()
        gene_info = defaultdict(list)
        for gi in qs.values('name', 'description', 'icon_css_class', gene_symbol_path).distinct():
            gene_symbol = gi[gene_symbol_path]
            gene_info[gene_symbol].append(gi)

        return Response(dict(gene_info))


class BatchGeneIdentifierForReleaseView(APIView):
    """ Needs to be a post as we can send a large number of genes """

    @extend_schema(
        summary="Map a batch of gene symbols to gene IDs for a GeneAnnotationRelease",
        request=OpenApiTypes.OBJECT,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, release_id):
        gene_annotation_release = get_object_or_404(GeneAnnotationRelease, pk=release_id)
        gene_symbols_json = request.data["gene_symbols_json"]
        gene_symbols = json.loads(gene_symbols_json)

        qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gene_annotation_release,
                                                  release_gene_symbol__gene_symbol__in=gene_symbols)
        gene_symbol_genes = defaultdict(list)
        for gene_symbol, gene_id in qs.values_list("release_gene_symbol__gene_symbol", "gene_id"):
            gene_symbol_genes[gene_symbol].append(gene_id)
        return Response(dict(gene_symbol_genes))


class SampleGeneListView(APIView):
    """ Update a SampleGeneList's active/visible status """

    @extend_schema(
        summary="Set active/visible status of a sample gene list",
        request=OpenApiTypes.OBJECT,
        responses=SampleGeneListSerializer,
    )
    def post(self, request, pk):
        sample_gene_list = get_object_or_404(SampleGeneList, pk=pk)
        sample_gene_list.sample.check_can_write(request.user)

        active = request.data.get("active")
        if active is not None:
            active = json.loads(active)
            if active:
                sample_gene_list.visible = True  # Active must be visible
                ActiveSampleGeneList.objects.update_or_create(sample=sample_gene_list.sample,
                                                              defaults={"sample_gene_list": sample_gene_list})

        visible = request.data.get("visible")
        if visible is not None:
            visible = json.loads(visible)
            # If you hide active gene list, it's no longer active...
            if visible is False:
                ActiveSampleGeneList.objects.filter(sample_gene_list=sample_gene_list).delete()

            sample_gene_list.visible = visible

        sample_gene_list.save()

        data = SampleGeneListSerializer(sample_gene_list, context={'request': request}).data
        return Response(data)


class _GenomeBuildContextMixin:
    """ Optional ?genome_build= restricts the gene versions returned to that build """

    def get_serializer_context(self):
        context = super().get_serializer_context()
        if genome_build_name := self.request.query_params.get("genome_build"):
            context["genome_build"] = GenomeBuild.get_name_or_alias(genome_build_name)
        return context


class GeneDetailView(_GenomeBuildContextMixin, RetrieveAPIView):
    """ A Gene, its per genome build versions, and the transcripts that place it on the genome """
    serializer_class = GeneDetailSerializer
    queryset = Gene.objects.all()
    lookup_url_kwarg = "gene_id"

    @extend_schema(
        summary="Retrieve a gene, its versions and transcripts",
        parameters=[OpenApiParameter(name="genome_build", type=OpenApiTypes.STR, location=OpenApiParameter.QUERY,
                                     description="Restrict versions to this genome build, eg 'GRCh38'")],
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)


class GeneSymbolDetailView(_GenomeBuildContextMixin, RetrieveAPIView):
    """ A gene symbol, its aliases, and the genes (RefSeq + Ensembl) assigned to it """
    serializer_class = GeneSymbolDetailSerializer
    queryset = GeneSymbol.objects.all()
    lookup_url_kwarg = "gene_symbol"

    @extend_schema(
        summary="Retrieve a gene symbol, its aliases, genes and transcripts",
        parameters=[OpenApiParameter(name="genome_build", type=OpenApiTypes.STR, location=OpenApiParameter.QUERY,
                                     description="Restrict versions to this genome build, eg 'GRCh38'")],
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)
