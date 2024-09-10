from rest_framework import serializers

from genes.serializers import GeneListSerializer
from seqauto.models import EnrichmentKit



class EnrichmentKitSerializer(serializers.ModelSerializer):
    gene_list = GeneListSerializer()
    enrichment_kit_type = serializers.SerializerMethodField()
    manufacturer = serializers.StringRelatedField()
    __str__ = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = EnrichmentKit
        fields = ('pk', 'name', 'version', 'enrichment_kit_type', 'manufacturer', 'gene_list', '__str__')

    def to_internal_value(self, data):
        # When POSTing, you can pass just name/version as a shortcut
        name = data.get('name')
        version = data.get('version')
        if enrichment_kit := EnrichmentKit.objects.filter(name=name, version=version).first():
            return enrichment_kit
        return super().to_internal_value(data)

    def get_enrichment_kit_type(self, obj):
        return obj.get_enrichment_kit_type_display()

    def get___str__(self, obj):
        return str(obj)


class EnrichmentKitSummarySerializer(serializers.ModelSerializer):
    """ Doesn't return genes (much faster)  """
    enrichment_kit_type = serializers.SerializerMethodField()
    manufacturer = serializers.StringRelatedField()
    __str__ = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = EnrichmentKit
        fields = ('pk', 'name', 'version', 'enrichment_kit_type', 'manufacturer', '__str__')

    def get_enrichment_kit_type(self, obj):
        return obj.get_enrichment_kit_type_display()

    def get___str__(self, obj):
        return str(obj)