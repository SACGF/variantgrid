from django.db import models
from rest_framework import serializers

from analysis.models import (
    AlleleFrequencyNode,
    AllVariantsNode,
    Analysis,
    AnalysisVariable,
    BuiltInFilterNode,
    CandidateSearchRun,
    ClassificationsNode,
    ClinVarNode,
    CohortNode,
    DamageNode,
    FilterNode,
    FilterNodeItem,
    GeneListNode,
    GeneListNodeGeneList,
    IntersectionNode,
    MergeNode,
    MOINode,
    MOINodeModeOfInheritance,
    MOINodeOntologyTerm,
    MOINodeSubmitter,
    PedigreeNode,
    PhenotypeNode,
    PhenotypeNodeOntologyTerm,
    PopulationNode,
    PopulationNodeGnomADPopulation,
    QuadNode,
    SampleNode,
    SelectedInParentNode,
    TagNode,
    TissueNode,
    TrioNode,
    VennNode,
    ZygosityNode,
)
from analysis.models.models_variant_tag import VariantTag
from analysis.models.nodes.analysis_node import (
    AnalysisNode,
    NodeAlleleFrequencyFilter,
    NodeAlleleFrequencyRange,
    NodeWiki,
)
from analysis.models.nodes.filters.conservation_node import ConservationNode
from genes.serializers import GeneListSerializer
from library.django_utils import get_model_fields
from library.django_utils.django_rest_utils import DynamicFieldsModelSerializer
from ontology.models import OntologyTerm
from ontology.serializers import OntologyTermSerializer
from snpdb.serializers import TimestampField, UserSerializer


class NodeAlleleFrequencyRangeSerializer(serializers.ModelSerializer):

    class Meta:
        model = NodeAlleleFrequencyRange
        fields = ('min', 'max')


class NodeAlleleFrequencyFilterSerializer(serializers.ModelSerializer):
    nodeallelefrequencyrange_set = NodeAlleleFrequencyRangeSerializer(many=True)

    class Meta:
        model = NodeAlleleFrequencyFilter
        fields = ('group_operation', 'nodeallelefrequencyrange_set')


class AnalysisVariableSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalysisVariable
        exclude = ("node", )  # So it works across systems and it's inside node anyway


class NodeWikiSerializer(serializers.ModelSerializer):
    class Meta:
        model = NodeWiki
        exclude = ("last_edited_by", )  # User PK - AnalysisNodeSerializer.create sets the importing user


class AnalysisSerializer(DynamicFieldsModelSerializer):
    class Meta:
        model = Analysis
        fields = "__all__"


def _analysis_node_fields(model):
    NODE_FIELDS = ["analysisvariable_set", "nodeallelefrequencyfilter", "nodewiki"]
    NODE_EXCLUDE = ["analysisnode_ptr", "children"]
    return get_model_fields(model, ignore_fields=NODE_EXCLUDE) + NODE_FIELDS


class AnalysisNodeSerializer(DynamicFieldsModelSerializer):
    analysisvariable_set = AnalysisVariableSerializer(many=True, allow_null=True)
    nodeallelefrequencyfilter = NodeAlleleFrequencyFilterSerializer(allow_null=True)
    # nodevcffilter_set = NodeVCFFilterSerializer(many=True)  # Don't do this as not consistent across VCFs
    nodewiki = NodeWikiSerializer(allow_null=True)

    class Meta:
        model = AnalysisNode
        fields = _analysis_node_fields(model)

    def create(self, validated_data):
        analysisvariable_set_data = validated_data.pop('analysisvariable_set')
        nodeallelefrequencyfilter_data = validated_data.pop('nodeallelefrequencyfilter', None)
        nodewiki_data = validated_data.pop('nodewiki')

        node = self.Meta.model.objects.create(**validated_data)

        for av_data in analysisvariable_set_data:
            AnalysisVariable.objects.create(node=node, **av_data)

        if nodeallelefrequencyfilter_data:
            nodeallelefrequencyrange_set_data = nodeallelefrequencyfilter_data.pop("nodeallelefrequencyrange_set")
            naf_filter, _ = NodeAlleleFrequencyFilter.objects.get_or_create(node=node, defaults=nodeallelefrequencyfilter_data)
            for filter_range_data in nodeallelefrequencyrange_set_data:
                NodeAlleleFrequencyRange.objects.create(filter=naf_filter, **filter_range_data)

        if nodewiki_data:
            NodeWiki.objects.create(node=node, **nodewiki_data)
        return node

    def get_extra_kwargs(self):
        """ A node the user hasn't configured yet holds Django's implicit '' for its unset text
            fields, and an export has to be able to bring that back """
        extra_kwargs = super().get_extra_kwargs()
        for field in self.Meta.model._meta.fields:
            if isinstance(field, (models.CharField, models.TextField)) and not field.blank:
                extra_kwargs.setdefault(field.name, {}).setdefault("allow_blank", True)
        return extra_kwargs

    @staticmethod
    def get_node_serializers() -> dict[str, 'AnalysisNodeSerializer']:
        node_serializers = {}
        for serializer_subclass in AnalysisNodeSerializer.__subclasses__():
            model_name = serializer_subclass.Meta.model._meta.label
            node_serializers[model_name] = serializer_subclass
        return node_serializers


class AlleleFrequencyNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = AlleleFrequencyNode
        fields = _analysis_node_fields(model)


class AllVariantsNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = AllVariantsNode
        fields = _analysis_node_fields(model)


class BuiltInFilterNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = BuiltInFilterNode
        fields = _analysis_node_fields(model)


class ClassificationsNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = ClassificationsNode
        fields = _analysis_node_fields(model)


class ClinVarNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = ClinVarNode
        fields = _analysis_node_fields(model)


class CohortNodeSerializer(AnalysisNodeSerializer):
    # Doesn't do this as it doesn't work between systems
    # cohortnodezygosityfilter_set = CohortNodeZygosityFilterSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = CohortNode
        fields = _analysis_node_fields(model)


class ConservationNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = ConservationNode
        fields = _analysis_node_fields(model)


class DamageNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = DamageNode
        fields = _analysis_node_fields(model)


class FilterNodeItemSerializer(serializers.ModelSerializer):
    class Meta:
        model = FilterNodeItem
        exclude = ("filter_node",)


class FilterNodeSerializer(AnalysisNodeSerializer):
    filternodeitem_set = FilterNodeItemSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = FilterNode
        fields = _analysis_node_fields(FilterNode) + ["filternodeitem_set"]

    def create(self, validated_data):
        filternodeitem_set_data = validated_data.pop('filternodeitem_set')
        node = super().create(validated_data)
        for fni_data in filternodeitem_set_data:
            FilterNodeItem.objects.create(filter_node=node, **fni_data)
        return node


class GeneListNodeGeneListSerializer(serializers.ModelSerializer):
    gene_list = GeneListSerializer()

    class Meta:
        model = GeneListNodeGeneList
        exclude = ("gene_list_node",)

    def to_internal_value(self, data):
        gene_list_data = data.get('gene_list')
        if isinstance(gene_list_data, dict):
            gl_serializer = GeneListSerializer(data=gene_list_data, context=self.context)
            gl_serializer.is_valid(raise_exception=True)
            return {'gene_list': gl_serializer.save()}
        return super().to_internal_value(data)


class GeneListNodeSerializer(AnalysisNodeSerializer):
    genelistnodegenelist_set = GeneListNodeGeneListSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = GeneListNode
        fields = _analysis_node_fields(model) + ["genelistnodegenelist_set"]

    def create(self, validated_data):
        genelistnodegenelist_set_data = validated_data.pop('genelistnodegenelist_set')
        node = super().create(validated_data)
        for gln_data in genelistnodegenelist_set_data:
            GeneListNodeGeneList.objects.create(gene_list_node=node, **gln_data)
        return node


class IntersectionNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = IntersectionNode
        fields = _analysis_node_fields(model)


class MergeNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = MergeNode
        fields = _analysis_node_fields(model)


class MOINodeOntologyTermSerializer(serializers.ModelSerializer):
    ontology_term = OntologyTermSerializer()

    class Meta:
        model = MOINodeOntologyTerm
        exclude = ("node",)

    def to_internal_value(self, data):
        ontology_term_data = data.get('ontology_term')
        if isinstance(ontology_term_data, dict):
            ontology_id = ontology_term_data.get('id')
            if ontology_id:
                return {'ontology_term': OntologyTerm.objects.get(pk=ontology_id)}
        return super().to_internal_value(data)


class MOINodeModeOfInheritanceSerializer(serializers.ModelSerializer):
    class Meta:
        model = MOINodeModeOfInheritance
        exclude = ("node",)


class MOINodeSubmitterSerializer(serializers.ModelSerializer):
    class Meta:
        model = MOINodeSubmitter
        exclude = ("node",)


class MOINodeSerializer(AnalysisNodeSerializer):
    moinodeontologyterm_set = MOINodeOntologyTermSerializer(many=True)
    moinodemodeofinheritance_set = MOINodeModeOfInheritanceSerializer(many=True)
    moinodesubmitter_set = MOINodeSubmitterSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = MOINode
        fields = _analysis_node_fields(model) + ["moinodeontologyterm_set", "moinodemodeofinheritance_set",
                                                "moinodesubmitter_set"]

    def create(self, validated_data):
        moinodeontologyterm_set_data = validated_data.pop('moinodeontologyterm_set')
        moinodemodeofinheritance_set_data = validated_data.pop('moinodemodeofinheritance_set')
        moinodesubmitter_set_data = validated_data.pop('moinodesubmitter_set')

        node = super().create(validated_data)

        for ontology_data in moinodeontologyterm_set_data:
            MOINodeOntologyTerm.objects.create(node=node, **ontology_data)
        for moi_data in moinodemodeofinheritance_set_data:
            MOINodeModeOfInheritance.objects.create(node=node, **moi_data)
        for submitter_data in moinodesubmitter_set_data:
            MOINodeSubmitter.objects.create(node=node, **submitter_data)

        return node


class PedigreeNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = PedigreeNode
        fields = _analysis_node_fields(model)


class PhenotypeNodeOntologyTermSerializer(serializers.ModelSerializer):
    ontology_term = OntologyTermSerializer()

    class Meta:
        model = PhenotypeNodeOntologyTerm
        exclude = ("phenotype_node",)

    def to_internal_value(self, data):
        ontology_term_data = data.get('ontology_term')
        if isinstance(ontology_term_data, dict):
            ontology_id = ontology_term_data.get('id')
            if ontology_id:
                return {'ontology_term': OntologyTerm.objects.get(pk=ontology_id)}
        return super().to_internal_value(data)


class PhenotypeNodeSerializer(AnalysisNodeSerializer):
    phenotypenodeontologyterm_set = PhenotypeNodeOntologyTermSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = PhenotypeNode
        fields = _analysis_node_fields(model) + ["phenotypenodeontologyterm_set"]

    def create(self, validated_data):
        phenotypenodeontologyterm_set = validated_data.pop('phenotypenodeontologyterm_set')

        node = super().create(validated_data)

        for ontology_data in phenotypenodeontologyterm_set:
            PhenotypeNodeOntologyTerm.objects.create(phenotype_node=node, **ontology_data)

        return node


class PopulationNodeGnomADPopulationSerializer(serializers.ModelSerializer):
    class Meta:
        model = PopulationNodeGnomADPopulation
        exclude = ("population_node", )  # So it works across systems and it's inside node anyway


class PopulationNodeSerializer(AnalysisNodeSerializer):
    populationnodegnomadpopulation_set = PopulationNodeGnomADPopulationSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = PopulationNode
        fields = _analysis_node_fields(model) + ["populationnodegnomadpopulation_set"]

    def create(self, validated_data):
        populationnodegnomadpopulation_set_data = validated_data.pop('populationnodegnomadpopulation_set')

        node = super().create(validated_data)

        for pn_gnomad_data in populationnodegnomadpopulation_set_data:
            PopulationNodeGnomADPopulation.objects.create(population_node=node, **pn_gnomad_data)

        return node


class QuadNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = QuadNode
        fields = _analysis_node_fields(model)


class SampleNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = SampleNode
        fields = _analysis_node_fields(model)


class SelectedInParentNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = SelectedInParentNode
        fields = _analysis_node_fields(model)


class TagNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = TagNode
        fields = _analysis_node_fields(model)


class TissueNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = TissueNode
        fields = _analysis_node_fields(model)


class TrioNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = TrioNode
        fields = _analysis_node_fields(model)


class VennNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = VennNode
        fields = _analysis_node_fields(model)


class ZygosityNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = ZygosityNode
        fields = _analysis_node_fields(model)


class VariantTagSerializer(serializers.ModelSerializer):
    analysis = AnalysisSerializer()
    user = UserSerializer()
    created = TimestampField()
    can_write = serializers.SerializerMethodField()

    class Meta:
        model = VariantTag
        fields = "__all__"

    def get_can_write(self, obj):
        user = self.context['request'].user
        return obj.can_write(user)


class CandidateSearchRunSerializer(serializers.ModelSerializer):
    class Meta:
        model = CandidateSearchRun
        fields = "__all__"
