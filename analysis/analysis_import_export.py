import json

from django.core.serializers import serialize

from django.contrib.auth.models import User

from analysis.models.nodes.node_utils import reload_analysis_nodes
from analysis.serializers import AnalysisNodeSerializer, AnalysisSerializer
from analysis.models import Analysis, GenomeBuild, AnnotationVersion, AnalysisEdge


def analysis_export_to_dict(analysis: Analysis) -> dict:
    node_serializers = AnalysisNodeSerializer.get_node_serializers()
    NODE_EXCLUDE = ["created", "modified", "version", "appearance_version", "ready", "valid", "count",
                    "celery_task", "db_pid"]

    analysis_serializer = AnalysisSerializer(analysis, exclude=["id", "created", "modified", "user",
                                                                "default_sort_by_column", "custom_columns_collection",
                                                                "visible", "annotation_version"])

    analysis_data = {
        "analysis": {"model": analysis._meta.label, "pk": analysis.pk, "fields": dict(analysis_serializer.data)}
    }

    nodes = []
    for node in analysis.analysisnode_set.all().select_subclasses():
        model_name = node._meta.label
        serializer = node_serializers.get(model_name, AnalysisNodeSerializer)
        fields = dict(serializer(node, exclude=NODE_EXCLUDE).data)
        nodes.append({"model": model_name, "pk": node.pk, "fields": dict(fields)})

    analysis_data["nodes"] = nodes
    analysis_data["edges"] = json.loads(serialize('json', AnalysisEdge.objects.filter(child__analysis=analysis)))
    return analysis_data


def analysis_export_to_file(analysis: Analysis, file):
    analysis_data = analysis_export_to_dict(analysis)
    json.dump(analysis_data, file)


def analysis_import(user: User, genome_build: GenomeBuild, filename,
                    annotation_version: AnnotationVersion = None) -> Analysis:
    if annotation_version is None:
        annotation_version = AnnotationVersion.latest(genome_build)

    with open(filename) as f:
        analysis_json = json.loads(f.read())

    node_serializers = {}
    for serializer_subclass in AnalysisNodeSerializer.__subclasses__():
        model_name = serializer_subclass.Meta.model._meta.label
        node_serializers[model_name] = serializer_subclass

    analysis_record = analysis_json["analysis"]
    analysis_kwargs = analysis_record["fields"]
    analysis_kwargs["user"] = user
    analysis_kwargs["genome_build"] = genome_build
    analysis_kwargs["annotation_version"] = annotation_version
    analysis = Analysis.objects.create(**analysis_kwargs)

    print(f"Creating analysis: {analysis.pk}")
    old_new_map = {}

    for node_record in analysis_json["nodes"]:
        model_name = node_record["model"]
        data = node_record["fields"]
        old_pk = data.pop("id")
        data["analysis"] = analysis.pk

        ############### TODO: Remove this
        data.pop("sample", None)
        # need to remove node from analysisvariable
        if analysisvariable_set := data.get("analysisvariable_set"):
            for av in analysisvariable_set:
                av.pop("node", None)  # Will provide this
        ############### END REMOVE

        serializer = node_serializers[model_name]
        s = serializer(data=data)
        if s.is_valid():
            try:
                node = s.save()
            except:
                print(f"Failed to save node: {old_pk}")
                print(data)
                raise

            old_new_map[old_pk] = node.pk
        else:
            print("Invalid:")
            print(s.errors)

    for edge_record in analysis_json["edges"]:
        edge_fields = edge_record["fields"]
        old_parent_id = edge_fields["parent"]
        old_child_id = edge_fields["child"]

        parent_id = old_new_map.get(old_parent_id)
        child_id = old_new_map.get(old_child_id)

        if parent_id and child_id:
            AnalysisEdge.objects.create(parent_id=parent_id, child_id=child_id)
        else:
            print(f"SKIPPED EDGE Parent: {old_parent_id} to Child: {old_child_id}")

    reload_analysis_nodes(analysis.pk)
    return analysis
