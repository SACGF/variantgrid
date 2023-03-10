import json
from collections import defaultdict

from django.http.response import Http404
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework import permissions
from rest_framework.generics import ListAPIView
from rest_framework.generics import RetrieveAPIView
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from genes.gene_matching import GeneSymbolMatcher, ReleaseGeneMatcher
from genes.models import GeneInfo, GeneList, GeneAnnotationRelease, \
    ReleaseGeneSymbolGene, PanelAppServer, SampleGeneList, ActiveSampleGeneList, create_fake_gene_list
from genes.panel_app import get_panel_app_panel_as_gene_list_json
from genes.panel_app import get_panel_app_results_by_gene_symbol_json
from genes.serializers import GeneInfoSerializer, GeneListGeneSymbolSerializer, GeneListSerializer, \
    GeneAnnotationReleaseSerializer, SampleGeneListSerializer
from library.constants import WEEK_SECS
from library.django_utils.django_rest_utils import MultipleFieldLookupMixin
from library.guardian_utils import DjangoPermission
from library.log_utils import get_traceback
from snpdb.models.models_enums import ImportStatus


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

    def get(self, request, *args, **kwargs):
        panel_app_id = self.kwargs['pk']
        data = get_panel_app_panel_as_gene_list_json(panel_app_id)
        return Response(data)


class GeneListView(APIView):

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
    permission_classes = (permissions.IsAuthenticated, WriteGeneListPermission)

    def get_gene_list_modifications(self, data):
        modifications = json.loads(data["modifications"])
        gene_additions = modifications.get('add')
        gene_deletions = modifications.get('delete')
        return gene_additions, gene_deletions

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
        except:
            gene_list.error_message = get_traceback()
            import_status = ImportStatus.ERROR

        gene_list.import_status = import_status
        gene_list.save()

        return Response({"pk": gene_list.pk})


class TextToGeneListView(APIView):
    """ Text to gene list (doesn't actually save). Used to eg check gene names """

    def post(self, request, *args, **kwargs):
        # Needed to be post as ~500 genes exceeded GET limit of ~4k
        name = self.request.data['name']
        gene_list_text = self.request.data['gene_list_text']

        gene_list = create_fake_gene_list(name=name, user=request.user)
        serializer = GeneListSerializer(gene_list, context={"request": request})
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


class GeneInfoView(ListAPIView):
    serializer_class = GeneInfoSerializer

    def get_queryset(self):
        gene_symbol = self.kwargs['gene_symbol']
        return GeneInfo.get_for_gene_symbol(gene_symbol)


class BatchGeneInfoView(APIView):
    """ Needs to be a post as we can send a large number of genes
        returns {gene_symbol : [gene_info_dict1, gene_info_dict2]} """

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

    def post(self, request, release_id):
        gene_annotation_release = get_object_or_404(GeneAnnotationRelease, pk=release_id)
        gene_symbols_json = request.data["gene_symbols_json"]
        gene_symbols = json.loads(gene_symbols_json)

        gm = ReleaseGeneMatcher(gene_annotation_release)
        gm.match_unmatched_symbols(gene_symbols)

        qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gene_annotation_release,
                                                  release_gene_symbol__gene_symbol__in=gene_symbols)
        gene_symbol_genes = defaultdict(list)
        for gene_symbol, gene_id in qs.values_list("release_gene_symbol__gene_symbol", "gene_id"):
            gene_symbol_genes[gene_symbol].append(gene_id)
        return Response(dict(gene_symbol_genes))


class SampleGeneListView(APIView):

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

        data = SampleGeneListSerializer(sample_gene_list).data
        return Response(data)
