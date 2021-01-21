from typing import Dict

from rest_framework import serializers

from analysis.models import AnalysisVariable, FilterNode, FilterNodeItem, PhenotypeNode, BuiltInFilterNode, \
    ClassificationsNode, DamageNode, ExpressionNode, VennNode, ZygosityNode, TrioNode, CohortNode, PedigreeNode, \
    TissueNode, SelectedInParentNode, SampleNode, PopulationNode, PopulationNodeGnomADPopulation, GeneListNode, \
    IntersectionNode, GeneListNodeGeneList, Analysis, AlleleFrequencyNode, AllVariantsNode, TagNode, MergeNode, \
    PhenotypeNodeOntologyTerm
from analysis.models.nodes.analysis_node import NodeAlleleFrequencyRange, NodeAlleleFrequencyFilter, AnalysisNode, \
    NodeWiki
from genes.serializers import GeneListSerializer
from library.django_utils import get_model_fields
from library.django_utils.django_rest_utils import DynamicFieldsModelSerializer
from ontology.serializers import OntologyTermSerializer


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
        fields = "__all__"


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

    @staticmethod
    def get_node_serializers() -> Dict[str, 'AnalysisNodeSerializer']:
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


class CohortNodeSerializer(AnalysisNodeSerializer):
    # Doesn't do this as it doesn't work between systems
    # cohortnodezygosityfilter_set = CohortNodeZygosityFilterSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = CohortNode
        fields = _analysis_node_fields(model)


class DamageNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = DamageNode
        fields = _analysis_node_fields(model)


class ExpressionNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = ExpressionNode
        fields = _analysis_node_fields(model)


class FilterNodeItemSerializer(serializers.ModelSerializer):
    class Meta:
        model = FilterNodeItem
        fields = "__all__"


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
        fields = "__all__"


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


class PedigreeNodeSerializer(AnalysisNodeSerializer):
    class Meta(AnalysisNodeSerializer.Meta):
        model = PedigreeNode
        fields = _analysis_node_fields(model)


class PhenotypeNodeOntologyTermSerializer(serializers.ModelSerializer):
    ontology_term = OntologyTermSerializer()

    class Meta:
        model = PhenotypeNodeOntologyTerm
        fields = "__all__"


class PhenotypeNodeSerializer(AnalysisNodeSerializer):
    phenotypenodeontology_set = PhenotypeNodeOntologyTermSerializer(many=True)

    class Meta(AnalysisNodeSerializer.Meta):
        model = PhenotypeNode
        fields = _analysis_node_fields(model) + ["phenotypenodehpo_set", "phenotypenodeomim_set"]

    def create(self, validated_data):
        phenotypenodeontology_set = validated_data.pop('phenotypenodeontology_set')

        node = super().create(validated_data)

        for ontology_data in phenotypenodeontology_set:
            print(f"ontology_data: {ontology_data}")
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
