from django.http.response import Http404
from django.utils.decorators import method_decorator
from django.views.decorators.cache import never_cache
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework.generics import RetrieveAPIView
from rest_framework.response import Response
from rest_framework.views import APIView

from annotation.models.models import ManualVariantEntryCollection, VariantAnnotationVersion
from annotation.serializers import (
    ManualVariantEntryCollectionSerializer,
    VariantAnnotationSerializer,
)
from snpdb.models import GenomeBuild, Variant


class VariantAnnotationView(APIView):
    """Retrieve the latest VariantAnnotation for a variant identified by coordinate string."""

    @extend_schema(
        summary="Retrieve latest variant annotation for a variant in a genome build",
        parameters=[
            OpenApiParameter("genome_build_name", OpenApiTypes.STR, OpenApiParameter.PATH,
                             description="Genome build name or alias, e.g. 'GRCh37' or 'GRCh38'"),
            OpenApiParameter("variant_string", OpenApiTypes.STR, OpenApiParameter.PATH,
                             description="Variant as 'chrom-position-ref-alt', e.g. '1-169519049-T-C'"),
        ],
        responses=VariantAnnotationSerializer,
    )
    def get(self, request, *args, **kwargs):
        genome_build_name = self.kwargs['genome_build_name']
        variant_string = self.kwargs['variant_string']
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        variant = Variant.get_from_string(variant_string, genome_build)

        if variant is None:
            raise Http404(variant_string)

        vav = VariantAnnotationVersion.latest(genome_build)
        va = variant.variantannotation_set.get(version=vav)
        serializer = VariantAnnotationSerializer(va)
        return Response(serializer.data)


@method_decorator(never_cache, name='get')
class ManualVariantEntryCollectionView(RetrieveAPIView):
    serializer_class = ManualVariantEntryCollectionSerializer
    lookup_field = 'pk'

    def get_queryset(self):
        # Doesn't need any security as they're just variants
        return ManualVariantEntryCollection.objects.all()
