from rest_framework import serializers

from ontology.models import OntologyTerm


class OntologyTermSerializer(serializers.ModelSerializer):
    ontology_service = serializers.SerializerMethodField()

    class Meta:
        model = OntologyTerm
        fields = ["id", "ontology_service", "name", "definition"]

    def get_ontology_service(self, obj):
        return obj.get_ontology_service_display()
