from rest_framework import serializers

from classification.models import Classification
from snpdb.serializers import VariantSerializer


# Is this used?
class ClassificationSerializer(serializers.ModelSerializer):
    variant = VariantSerializer()
    clinical_significance = serializers.SerializerMethodField()
    lab = serializers.StringRelatedField()

    class Meta:
        model = Classification
        fields = '__all__'

    def get_clinical_significance(self, obj):
        return obj.get_clinical_significance_display()
