import numpy as np
from rest_framework import serializers

from genes.serializers import TranscriptSerializer, GeneSymbolSerializer, \
    TranscriptVersionSerializer
from seqauto.models import GoldCoverageSummary, GoldReference
from seqauto.serializers.enrichment_kit_serializers import EnrichmentKitSummarySerializer


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
