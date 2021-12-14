from rest_framework import serializers

from annotation.models.models import VariantAnnotationVersion, VariantAnnotation, DiseaseValidity, \
    ManualVariantEntryCollection
from ontology.serializers import OntologyTermSerializer
from snpdb.serializers import VariantSerializer


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

    class Meta:
        model = ManualVariantEntryCollection
        fields = ('created', 'import_status', 'first_entry_text', 'first_variant_id', 'first_variant_annotation_status')

    def get_first_entry_text(self, obj):
        entry_text = None
        if mve := obj.manualvariantentry_set.all().order_by("pk").first():
            entry_text = mve.entry_text
        return entry_text

    def get_first_variant_id(self, obj):
        variant_id = None
        if mve := obj.manualvariantentry_set.all().order_by("pk").first():
            if cve := mve.unique_created_variants.filter(variant__isnull=False).order_by("variant_id").first():
                variant_id = cve.variant_id
        return variant_id

    def get_first_variant_annotation_status(self, obj):
        return "unknown"
