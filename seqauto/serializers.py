

from rest_framework import serializers

from genes.serializers import GeneListSerializer, GeneSerializer, TranscriptSerializer, GeneSymbolSerializer
from seqauto.models import GoldCoverageSummary, GoldReference, EnrichmentKit
import numpy as np


class EnrichmentKitSerializer(serializers.ModelSerializer):
    gene_list = GeneListSerializer()
    enrichment_kit_type = serializers.SerializerMethodField()
    manufacturer = serializers.StringRelatedField()
    __str__ = serializers.SerializerMethodField()

    class Meta:
        model = EnrichmentKit
        fields = ('pk', 'name', 'version', 'enrichment_kit_type', 'manufacturer', 'gene_list', '__str__')

    def get_enrichment_kit_type(self, obj):
        return obj.get_enrichment_kit_type_display()

    def get___str__(self, obj):
        return str(obj)


class EnrichmentKitSummarySerializer(serializers.ModelSerializer):
    """ Doesn't return genes (much faster)  """
    enrichment_kit_type = serializers.SerializerMethodField()
    manufacturer = serializers.StringRelatedField()
    __str__ = serializers.SerializerMethodField()

    class Meta:
        model = EnrichmentKit
        fields = ('pk', 'name', 'version', 'enrichment_kit_type', 'manufacturer', '__str__')

    def get_enrichment_kit_type(self, obj):
        return obj.get_enrichment_kit_type_display()

    def get___str__(self, obj):
        return str(obj)


class GoldReferenceSerializer(serializers.ModelSerializer):
    enrichment_kit = EnrichmentKitSummarySerializer()

    class Meta:
        model = GoldReference
        fields = ('enrichment_kit', 'created')


class GoldCoverageSummarySerializer(serializers.ModelSerializer):
    gold_reference = GoldReferenceSerializer()
    gene_symbol = GeneSymbolSerializer()
    transcript = TranscriptSerializer()
    transcript_version = TranscriptVersionSerializer()
    standard_error = serializers.SerializerMethodField()

    class Meta:
        model = GoldCoverageSummary
        fields = '__all__'

    def get_standard_error(self, obj):
        """ This can occasionally be NaN which isn't valid JSON """
        standard_error = obj.standard_error
        if np.isnan(standard_error):
            standard_error = -1
        return standard_error
