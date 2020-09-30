from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from rest_framework import permissions, viewsets
from rest_framework.generics import RetrieveAPIView, get_object_or_404
from rest_framework.response import Response
from rest_framework.views import APIView

from library.constants import MINUTE_SECS
from patients.models_enums import Zygosity
from snpdb.clingen_allele import get_variant_allele_for_variant
from snpdb.models import Sample, Variant, Trio
from snpdb.models.models_vcf import Project
from snpdb.serializers import TrioSerializer, VariantAlleleSerializer, ProjectSerializer


class VariantZygosityForSampleView(APIView):
    """ Returns Zygosity for variant/sample ("." if missing) """

    def get(self, request, *args, **kwargs):
        variant_id = self.kwargs['variant_id']
        sample_id = self.kwargs['sample_id']

        variant = Variant.objects.get(pk=variant_id)
        sample = Sample.get_for_user(self.request.user, sample_id)

        if sample_genotype := sample.get_genotype(variant):
            zygosity = sample_genotype.zygosity
        else:
            zygosity = Zygosity.MISSING

        zygosity_dict = {"zygosity": Zygosity.display(zygosity)}
        return Response(zygosity_dict)


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class TrioView(RetrieveAPIView):
    serializer_class = TrioSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    lookup_fields = ('pk',)

    def get_queryset(self):
        return Trio.filter_for_user(self.request.user)


class VariantAlleleForVariantView(APIView):

    def get(self, request, *args, **kwargs):
        variant_id = self.kwargs['variant_id']
        variant = get_object_or_404(Variant, pk=variant_id)
        genome_build = variant.genome_build

        variant_allele = get_variant_allele_for_variant(genome_build, variant)
        data = VariantAlleleSerializer.data_with_link_data(variant_allele)
        return Response(data)


class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = ProjectSerializer
