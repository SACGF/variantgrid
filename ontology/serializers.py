from rest_framework import serializers

from ontology.models import OntologyTerm, OntologyTermRelation


class OntologyTermSerializer(serializers.ModelSerializer):
    ontology_service = serializers.SerializerMethodField()
    absolute_url = serializers.SerializerMethodField()

    class Meta:
        model = OntologyTerm
        fields = ["id", "ontology_service", "name", "definition", "absolute_url"]

    def get_ontology_service(self, obj: OntologyTerm):
        return obj.get_ontology_service_display()

    def get_absolute_url(self, obj: OntologyTerm):
        return obj.get_absolute_url()



class OntologyTermRelationSerializer(serializers.ModelSerializer):
    source_term = OntologyTermSerializer()
    dest_term = OntologyTermSerializer()

    class Meta:
        model = OntologyTermRelation
        fields = ["source_term", "dest_term", "extra"]
