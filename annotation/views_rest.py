from django.http.response import Http404
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework.generics import ListAPIView, RetrieveAPIView
from rest_framework.response import Response
from rest_framework.views import APIView

from annotation.models.models import DiseaseValidity, VariantAnnotationVersion, ManualVariantEntryCollection
from annotation.serializers import VariantAnnotationSerializer, DiseaseValiditySerializer, \
    ManualVariantEntryCollectionSerializer
from library.constants import WEEK_SECS
from snpdb.models import GenomeBuild, Variant


@method_decorator(cache_page(WEEK_SECS), name='get')
class GeneDiseaseValidityView(ListAPIView):
    serializer_class = DiseaseValiditySerializer

    def get_queryset(self):
        gene_symbol = self.kwargs['gene_symbol']
        return DiseaseValidity.objects.filter(genediseasevalidity__gene_symbol=gene_symbol).distinct()


class VariantAnnotationView(APIView):

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


class ManualVariantEntryCollectionView(RetrieveAPIView):
    serializer_class = ManualVariantEntryCollectionSerializer
    lookup_field = 'pk'

    def get_queryset(self):
        # Doesn't need any security as they're just variants
        return ManualVariantEntryCollection.objects.all()
