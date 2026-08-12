"""
#1707 - the patient / specimen / extraction API, so a lab client can accession before or alongside
posting a run.

Every serializer here is an upsert keyed on the identifiers the client sent (@see
patients.external_references), so a client re-posting a run gets the same rows back rather than
duplicates.
"""
from typing import Optional

from django.conf import settings
from django.db import transaction
from django.db.models import Model, Q
from drf_spectacular.utils import extend_schema_field
from rest_framework import serializers

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.external_references import (
    ExternalReference,
    ExternalReferenceError,
    resolve_reference,
)
from patients.models import (
    Extraction,
    ExternalModelManager,
    ExternalPK,
    Patient,
    Specimen,
    SpecimenMeasure,
)


def reference_json(obj: Optional[Model]) -> Optional[dict]:
    """ The identifiers a client could name an object by - what the API echoes back for a relation """
    if obj is None:
        return None
    data = {}
    if local_reference := getattr(obj, obj.LOCAL_REFERENCE_FIELD, None):
        data["reference_id"] = local_reference
    if external_pk := obj.external_pk:
        data.update(code=external_pk.code, external_type=external_pk.external_type,
                    external_manager=external_pk.external_manager_id)
    return data


@extend_schema_field({
    "oneOf": [
        {"type": "string", "description": "The local reference (reference_id / patient_code)"},
        {"type": "object",
         "properties": {
             "reference_id": {"type": "string"},
             "code": {"type": "string", "description": "ExternalPK.code - send with external_type"},
             "external_type": {"type": "string"},
             "external_manager": {"type": "string"},
         }},
    ],
})
class ExternalReferenceField(serializers.Field):
    """ The bare-string or object form of an identifier - @see patients.external_references """

    def to_internal_value(self, data) -> ExternalReference:
        try:
            return ExternalReference.from_data(data)
        except ExternalReferenceError as e:
            raise serializers.ValidationError(str(e)) from e

    def to_representation(self, value):
        if isinstance(value, ExternalReference):
            return value.as_json()
        if isinstance(value, Model):
            return reference_json(value)
        return value


def resolve_or_raise(model: type[Model], reference: ExternalReference, user) -> Model:
    """ A row that has nowhere to live without its parent gets a 400 rather than a parked reference -
        unlike a VCF naming an extraction, which has a row to park on """
    resolved = resolve_reference(model, reference, user)
    if not resolved.matched:
        raise serializers.ValidationError(resolved.error)
    return resolved.obj


class ExternalPKSerializer(serializers.Serializer):
    """ All three of code / external_type / external_manager are required, because that triple is what
        ExternalPK is unique on - anything less names a row we could not look up again. """
    code = serializers.CharField()
    external_type = serializers.CharField()
    external_manager = serializers.CharField()

    def validate_external_manager(self, value) -> ExternalModelManager:
        """ On an intranet deployment the set of managers is known, so an unfamiliar name is a typo
            rather than a new system - and silently creating one would decide can_modify (ie whether
            VG users may edit those records) as a side effect of a misspelling. A superuser creating
            one through the API is a deliberate act, so that is allowed either way """
        if manager := ExternalModelManager.objects.filter(name=value).first():
            return manager

        user = self.context["request"].user
        if settings.PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY and not user.is_superuser:
            known = ", ".join(ExternalModelManager.objects.values_list("name", flat=True)) or "(none configured)"
            raise serializers.ValidationError(f"Unknown external_manager '{value}'. Known: {known}")
        return ExternalModelManager.objects.create(name=value)

    @staticmethod
    def get_or_create(validated_data) -> ExternalPK:
        external_pk, _ = ExternalPK.objects.get_or_create(**validated_data)
        return external_pk

    def to_representation(self, instance) -> Optional[dict]:
        if instance is None:
            return None
        return {"code": instance.code, "external_type": instance.external_type,
                "external_manager": instance.external_manager_id}


class ExternallyManagedModelSerializer(serializers.ModelSerializer):
    """ Upsert keyed on the identifiers, so a client that re-posts a run gets the same rows """
    external_pk = ExternalPKSerializer(required=False)

    def _local_reference_q(self, validated_data) -> Optional[Q]:
        """ How a re-post finds the row it made last time, where no ExternalPK was sent """
        raise NotImplementedError

    def _post_create(self, instance, user):
        """ Anything a brand new row needs beyond being saved """

    def _get_existing(self, validated_data, external_pk_data, user) -> Optional[Model]:
        qs = self.Meta.model.filter_for_user(user)
        if external_pk_data:
            existing = qs.filter(external_pk__code=external_pk_data["code"],
                                 external_pk__external_type=external_pk_data["external_type"],
                                 external_pk__external_manager=external_pk_data["external_manager"]).first()
            if existing:
                return existing
        if q := self._local_reference_q(validated_data):
            return qs.filter(q).first()
        return None

    def create(self, validated_data):
        external_pk_data = validated_data.pop("external_pk", None)
        user = self.context["request"].user
        with transaction.atomic():
            instance = self._get_existing(validated_data, external_pk_data, user)
            created = instance is None
            if created:
                instance = self.Meta.model(**validated_data)
            else:
                instance.check_can_write(user)
                for field, value in validated_data.items():
                    setattr(instance, field, value)
            if external_pk_data:
                instance.external_pk = ExternalPKSerializer.get_or_create(external_pk_data)
            instance.save()
            if created:
                self._post_create(instance, user)
        return instance

    def update(self, instance, validated_data):
        external_pk_data = validated_data.pop("external_pk", None)
        user = self.context["request"].user
        instance.check_can_write(user)
        if external_pk_data:
            instance.external_pk = ExternalPKSerializer.get_or_create(external_pk_data)
        return super().update(instance, validated_data)


class PatientSerializer(ExternallyManagedModelSerializer):
    class Meta:
        model = Patient
        fields = ("id", "patient_code", "family_code", "first_name", "last_name", "date_of_birth",
                  "date_of_death", "sex", "affected", "external_pk")
        read_only_fields = ("id",)

    def _local_reference_q(self, validated_data) -> Optional[Q]:
        if patient_code := validated_data.get("patient_code"):
            return Q(patient_code=patient_code)
        return None

    def _post_create(self, instance, user):
        assign_permission_to_user_and_groups(user, instance)


class SpecimenSerializer(ExternallyManagedModelSerializer):
    patient = ExternalReferenceField()

    class Meta:
        model = Specimen
        fields = ("id", "patient", "reference_id", "description", "collected_by", "collection_date",
                  "received_date", "tissue_status", "external_pk")
        read_only_fields = ("id",)

    def validate_patient(self, reference: ExternalReference) -> Patient:
        """ A specimen has nowhere to live without its patient, so an unresolvable one is a 400 rather
            than a parked row - unlike a VCF naming an extraction, which has a row to park on """
        return resolve_or_raise(Patient, reference, self.context["request"].user)

    def _local_reference_q(self, validated_data) -> Optional[Q]:
        # reference_id is unique per patient rather than globally, so both name the row
        return Q(patient=validated_data["patient"], reference_id=validated_data["reference_id"])


class ExtractionSerializer(ExternallyManagedModelSerializer):
    specimen = ExternalReferenceField()

    class Meta:
        model = Extraction
        fields = ("id", "specimen", "reference_id", "nucleic_acid_source", "extraction_date",
                  "external_pk")
        read_only_fields = ("id",)

    def validate_specimen(self, reference: ExternalReference) -> Specimen:
        return resolve_or_raise(Specimen, reference, self.context["request"].user)

    def _local_reference_q(self, validated_data) -> Optional[Q]:
        if reference_id := validated_data.get("reference_id"):
            return Q(specimen=validated_data["specimen"], reference_id=reference_id)
        return None


class SpecimenMeasureFieldsSerializer(serializers.ModelSerializer):
    """ A measure without its specimen - the bulk call names that once, for the whole list """
    extraction = ExternalReferenceField(required=False)

    class Meta:
        model = SpecimenMeasure
        fields = ("id", "extraction", "measure_type", "value", "unit", "call", "threshold",
                  "threshold_source", "method", "source_payload", "measured_date")
        read_only_fields = ("id",)

    def validate_extraction(self, reference: ExternalReference) -> Extraction:
        return resolve_or_raise(Extraction, reference, self.context["request"].user)


class SpecimenMeasureSerializer(SpecimenMeasureFieldsSerializer):
    """ #1559 - a scalar measured on the material, transcribed by the client out of vendor output """
    specimen = ExternalReferenceField()

    class Meta(SpecimenMeasureFieldsSerializer.Meta):
        fields = ("specimen",) + SpecimenMeasureFieldsSerializer.Meta.fields

    def validate_specimen(self, reference: ExternalReference) -> Specimen:
        return resolve_or_raise(Specimen, reference, self.context["request"].user)

    def create(self, validated_data):
        """ A measure has no row to live on without its specimen, so an unresolvable one is a 400 -
            unlike a VCF naming an extraction, which parks on the sample it already created.

            A re-post replaces: there is one current MSI for a specimen, and the report wants that
            rather than a history to choose from. TimeStampedModel.modified records when it changed """
        return upsert_specimen_measure(validated_data.pop("specimen"), validated_data,
                                       self.context["request"].user)

    def update(self, instance, validated_data):
        instance.specimen.check_can_write(self.context["request"].user)
        validated_data["user"] = self.context["request"].user
        return super().update(instance, validated_data)


def upsert_specimen_measure(specimen: Specimen, validated_data: dict, user) -> SpecimenMeasure:
    specimen.check_can_write(user)
    defaults = {**validated_data, "specimen": specimen, "user": user}
    measure, _ = SpecimenMeasure.objects.update_or_create(specimen=specimen,
                                                          measure_type=validated_data["measure_type"],
                                                          defaults=defaults)
    return measure


class SpecimenMeasureBulkCreateSerializer(serializers.Serializer):
    """ A run's TMB, MSI, GIS, tumour fraction and ploidy post in one call, all against one specimen """
    specimen = ExternalReferenceField()
    measures = SpecimenMeasureFieldsSerializer(many=True)

    def validate_specimen(self, reference: ExternalReference) -> Specimen:
        return resolve_or_raise(Specimen, reference, self.context["request"].user)

    def create(self, validated_data):
        specimen = validated_data["specimen"]
        user = self.context["request"].user
        measures = []
        with transaction.atomic():
            for measure_data in validated_data["measures"]:
                measures.append(upsert_specimen_measure(specimen, measure_data, user))
        return {"specimen": specimen, "measures": measures}
