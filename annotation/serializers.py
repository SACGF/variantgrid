from rest_framework import serializers

from annotation.models import AnnotationStatus
from annotation.models.models import VariantAnnotationVersion, VariantAnnotation, DiseaseValidity, \
    ManualVariantEntryCollection
from ontology.serializers import OntologyTermSerializer
from snpdb.serializers import VariantSerializer, TimestampField


class DiseaseValiditySerializer(serializers.ModelSerializer):
    ontology_term = OntologyTermSerializer()
    classification = serializers.SerializerMethodField()
    gene_disease_curator = serializers.StringRelatedField()

    class Meta:
        model = DiseaseValidity
        fields = '__all__'

    def get_classification(self, obj):
        return obj.get_classification_display()


class VariantAnnotationVersionSerializer(serializers.ModelSerializer):

    class Meta:
        model = VariantAnnotationVersion
        fields = ('pk', 'annotation_date')


class VariantAnnotationSerializer(serializers.ModelSerializer):
    version = VariantAnnotationVersionSerializer()
    variant = VariantSerializer()

    class Meta:
        model = VariantAnnotation
        exclude = ('annotation_run',)


class ManualVariantEntryCollectionSerializer(serializers.ModelSerializer):
    first_entry_text = serializers.SerializerMethodField()
    first_variant_id = serializers.SerializerMethodField()
    first_variant_annotation_status = serializers.SerializerMethodField()
    created = TimestampField()
    is_ready = serializers.SerializerMethodField()

    class Meta:
        model = ManualVariantEntryCollection
        fields = '__all__'

    def get_first_entry_text(self, obj: ManualVariantEntryCollection):
        entry_text = None
        if mve := obj.first_entry:
            entry_text = mve.entry_text
        return entry_text

    def get_first_variant_id(self, obj: ManualVariantEntryCollection):
        variant_id = None
        if variant := obj.first_variant:
            variant_id = variant.pk
        return variant_id

    def get_first_variant_annotation_status(self, obj: ManualVariantEntryCollection):
        annotation_status = "Creating Variant"
        if obj.first_variant:
            annotation_status = "Variant Created"

        if ar := obj.first_variant_annotation_run():
            annotation_status = ar.get_status_display()
        return annotation_status

    def get_is_ready(self, obj: ManualVariantEntryCollection):
        ready = False
        if variant := obj.first_variant:
            if variant.is_reference:
                ready = True  # Won't be annotated

        if ar := obj.first_variant_annotation_run():
            ready = ar.status == AnnotationStatus.FINISHED
        return ready
