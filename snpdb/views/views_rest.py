from django.conf import settings
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema, OpenApiParameter
from rest_framework import permissions, viewsets
from rest_framework.exceptions import ValidationError
from rest_framework.generics import RetrieveAPIView, get_object_or_404
from rest_framework.response import Response
from rest_framework.views import APIView

from library.constants import MINUTE_SECS
from patients.models_enums import Zygosity
from snpdb.clingen_allele import get_variant_allele_for_variant
from snpdb.models import Sample, Variant, Trio, Quad, GenomeBuild
from snpdb.models.models_vcf import Project
from snpdb.serializers import QuadSerializer, TrioSerializer, VariantAlleleSerializer, ProjectSerializer
from snpdb.variant_sample_information import VariantSampleGenotypes


class VariantZygosityForSampleView(APIView):
    """ Returns Zygosity for variant/sample ("." if missing) """

    @extend_schema(
        summary="Get zygosity for a variant in a sample",
        responses=OpenApiTypes.OBJECT,
    )
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


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class QuadView(RetrieveAPIView):
    serializer_class = QuadSerializer
    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)
    lookup_fields = ('pk',)

    def get_queryset(self):
        return Quad.filter_for_user(self.request.user)


class VariantAlleleForVariantView(APIView):
    """ Returns the VariantAllele for a variant in a genome build, with variant link data """

    @extend_schema(
        summary="Get VariantAllele (with link data) for a variant in a genome build",
        responses=VariantAlleleSerializer,
        parameters=[
            OpenApiParameter("genome_build_name", OpenApiTypes.STR, OpenApiParameter.PATH,
                             description='Genome build name or alias, e.g. "GRCh37" or "GRCh38"'),
        ],
    )
    def get(self, request, *args, **kwargs):
        variant = get_object_or_404(Variant, pk=self.kwargs['variant_id'])
        genome_build = GenomeBuild.get_name_or_alias(self.kwargs['genome_build_name'])

        variant_allele = get_variant_allele_for_variant(genome_build, variant)
        data = VariantAlleleSerializer.data_with_link_data(variant_allele)
        return Response(data)


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class VariantSampleGenotypesView(APIView):
    """ Samples/genotypes for a variant, scoped to the samples the user can read.

        Keyed on variant (not variant + build) - a variant knows its builds via its contig, and
        builds sharing a contig (MT, unplaced scaffolds) resolve to the same variant. Each row is
        labelled with its own sample's VCF build.

        Cached briefly and keyed on the session (permission scoped payload) - a longer TTL would hide
        newly imported VCFs and newly published classifications from someone refreshing to see them. """

    @extend_schema(
        summary="Get samples/genotypes for a variant",
        responses=OpenApiTypes.OBJECT,
        parameters=[
            OpenApiParameter("limit", OpenApiTypes.INT, OpenApiParameter.QUERY,
                             description="Maximum rows returned (0 for all). Counts are always exact. "
                                         "Defaults to settings.VARIANT_DETAILS_SAMPLES_MAX_ROWS"),
        ],
    )
    def get(self, request, *args, **kwargs):
        variant = get_object_or_404(Variant, pk=self.kwargs['variant_id'])
        # limit=0 is "give me everything", which is max_rows=None (0 there means materialise nothing)
        max_rows = self._get_limit(request) or None
        # Stop scanning once the grid is full - "load all" (limit=0) is what does the exact count
        variant_sample_genotypes = VariantSampleGenotypes(request.user, variant, max_rows=max_rows,
                                                          stop_when_full=max_rows is not None)
        return Response(variant_sample_genotypes.to_json())

    @staticmethod
    def _get_limit(request) -> int:
        if (limit_param := request.query_params.get("limit")) is None:
            return settings.VARIANT_DETAILS_SAMPLES_MAX_ROWS
        try:
            limit = int(limit_param)
        except ValueError as ve:
            raise ValidationError({"limit": f"'{limit_param}' is not an integer"}) from ve
        if limit < 0:
            raise ValidationError({"limit": "must be 0 (all rows) or greater"})
        return limit


class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = ProjectSerializer
