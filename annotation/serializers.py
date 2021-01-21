from rest_framework import serializers

from annotation.models.models import EnsemblGeneAnnotationVersion, \
    EnsemblGeneAnnotation, VariantAnnotationVersion, VariantAnnotation, \
    DiseaseValidity
from annotation.models.models_mim_hpo import HumanPhenotypeOntology, HPOSynonym, \
    MIMMorbid, MIMMorbidAlias
from genes.serializers import GeneSerializer
from ontology.serializers import OntologyTermSerializer
from snpdb.serializers import VariantSerializer


class EnsemblGeneAnnotationVersionSerializer(serializers.ModelSerializer):

    class Meta:
        model = EnsemblGeneAnnotationVersion
        fields = ('pk', 'annotation_date')


class EnsemblGeneAnnotationSerializer(serializers.ModelSerializer):
    version = EnsemblGeneAnnotationVersionSerializer()
    gene = GeneSerializer()

    class Meta:
        model = EnsemblGeneAnnotation
        fields = '__all__'


class HumanPhenotypeOntologySerializer(serializers.ModelSerializer):

    class Meta:
        model = HumanPhenotypeOntology
        fields = ("name", "definition")


class HPOSynonymSerializer(serializers.ModelSerializer):
    hpo = HumanPhenotypeOntologySerializer()

    class Meta:
        model = HPOSynonym
        fields = '__all__'


class MIMMorbidSerializer(serializers.ModelSerializer):

    class Meta:
        model = MIMMorbid
        fields = '__all__'


class MIMMorbidAliasSerializer(serializers.ModelSerializer):
    mim_morbid = MIMMorbidSerializer()

    class Meta:
        model = MIMMorbidAlias
        fields = '__all__'


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
