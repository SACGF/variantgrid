from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import OpenApiParameter, extend_schema
from rest_framework.exceptions import NotFound
from rest_framework.generics import RetrieveAPIView, get_object_or_404
from rest_framework.response import Response
from rest_framework.views import APIView

from pathtests.models import PathologyTest, PathologyTestVersion
from pathtests.serializers import PathologyTestVersionSerializer


class PathologyTestVersionView(RetrieveAPIView):
    serializer_class = PathologyTestVersionSerializer
    lookup_field = 'pk'

    def get_queryset(self):
        return PathologyTestVersion.objects.all()


class PathologyTestLatestVersionView(APIView):
    """ Retrieves the last *confirmed* pathology test version for a pathology test """

    #@method_decorator(cache_page(WEEK_SECS))
    @extend_schema(
        summary="Retrieve the latest confirmed version of a pathology test by name",
        parameters=[
            OpenApiParameter("name", OpenApiTypes.STR, OpenApiParameter.PATH,
                             description="Pathology test name"),
        ],
        responses=PathologyTestVersionSerializer,
    )
    def get(self, request, *args, **kwargs):
        name = self.kwargs['name']

        pathology_test = get_object_or_404(PathologyTest, name=name)
        ptv = pathology_test.get_active_test_version()
        if not ptv:
            raise NotFound(f"No confirmed pathology test for {name}")

        serializer = PathologyTestVersionSerializer(ptv, context={"request": request})
        return Response(serializer.data)
