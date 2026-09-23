from rest_framework import viewsets

from patients.models import Extraction, Patient, Specimen
from patients.serializers import (
    ExtractionSerializer,
    PatientSerializer,
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
