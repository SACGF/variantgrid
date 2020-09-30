from rest_framework import serializers

from snpdb.serializers import VariantSerializer
from variantclassification.models import VariantClassification


class VariantClassificationSerializer(serializers.ModelSerializer):
    variant = VariantSerializer()
    clinical_significance = serializers.SerializerMethodField()
    lab = serializers.StringRelatedField()

    class Meta:
        model = VariantClassification
        fields = '__all__'

    def get_clinical_significance(self, obj):
        return obj.get_clinical_significance_display()
