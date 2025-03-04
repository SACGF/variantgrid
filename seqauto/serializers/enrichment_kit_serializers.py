from rest_framework import serializers

from genes.serializers import GeneListSerializer
from seqauto.models import EnrichmentKit



class EnrichmentKitSerializer(serializers.ModelSerializer):
    gene_list = GeneListSerializer(read_only=True)
    enrichment_kit_type = serializers.SerializerMethodField(read_only=True)
    manufacturer = serializers.StringRelatedField(read_only=True)
    __str__ = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = EnrichmentKit
        fields = ('pk', 'name', 'version', 'enrichment_kit_type', 'manufacturer', 'gene_list', '__str__')

    def to_internal_value(self, data):
        # When POSTing, you can pass just name/version as a shortcut
        name = data.get('name')
        version = data.get('version')
        if enrichment_kit := EnrichmentKit.objects.filter(name=name, version=version).first():
            data = self.to_representation(enrichment_kit)
        return super().to_internal_value(data)

    def create(self, validated_data):
        name = validated_data.get('name')
        version = validated_data.get('version')
        instance, _created = EnrichmentKit.objects.get_or_create(
            name=name,
            version=version,
            defaults=validated_data
        )
        return instance

    @staticmethod
    def get_from_data(data):
        return EnrichmentKit.objects.get(name=data["name"], version=data["version"])

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