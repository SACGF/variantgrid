from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema
from rest_framework import status, viewsets
from rest_framework.response import Response
from rest_framework.views import APIView

from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.serializers import (
    ExtractionSerializer,
    PatientSerializer,
    SpecimenMeasureBulkCreateSerializer,
    SpecimenMeasureSerializer,
    SpecimenSerializer,
)
from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions


class PatientViewSet(viewsets.ModelViewSet):
    serializer_class = PatientSerializer

    def get_queryset(self):
        return Patient.filter_for_user(self.request.user)


class SpecimenViewSet(viewsets.ModelViewSet):
    serializer_class = SpecimenSerializer

    def get_queryset(self):
        return Specimen.filter_for_user(self.request.user)


class ExtractionViewSet(viewsets.ModelViewSet):
    serializer_class = ExtractionSerializer

    def get_queryset(self):
        return Extraction.filter_for_user(self.request.user)

    def perform_create(self, serializer):
        super().perform_create(serializer)
        # New upstream data - re-run anything parked waiting for exactly this
        reconcile_pending_extractions.delay()


class SpecimenMeasureViewSet(viewsets.ModelViewSet):
    serializer_class = SpecimenMeasureSerializer

    def get_queryset(self):
        return SpecimenMeasure.filter_for_user(self.request.user)


class SpecimenMeasureBulkCreateView(APIView):
    """ A run's measures (TMB, MSI, GIS, tumour fraction, ploidy) against one specimen in one call """

    @extend_schema(
        summary="Bulk create specimen measures for one specimen",
        request=SpecimenMeasureBulkCreateSerializer,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, *args, **kwargs):
        serializer = SpecimenMeasureBulkCreateSerializer(data=request.data,
                                                         context={"request": request})
        serializer.is_valid(raise_exception=True)
        result = serializer.save()
        response = {
            "specimen": str(result["specimen"]),
            "measures": [str(measure) for measure in result["measures"]],
        }
        return Response(response, status=status.HTTP_201_CREATED)
