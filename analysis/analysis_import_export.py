import json
import logging
from typing import Optional

from django.contrib.auth.models import User
from django.core.serializers import serialize
from django.db import transaction

from analysis.enums import AnalysisTemplateType
from analysis.models import (
    Analysis,
    AnalysisEdge,
    AnalysisNodeCountConfiguration,
    AnnotationVersion,
    GenomeBuild,
    VennNode,
)
from analysis.models.nodes.node_utils import reload_analysis_nodes
from analysis.serializers import AnalysisNodeSerializer, AnalysisSerializer
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.enums import ZygosityRequirement

# References into a server's own data - an exported analysis arrives with its source nodes
# unconfigured, like a template before its variables are populated
NON_PORTABLE_NODE_FIELDS = [
    "cloned_from",
    "cohort",
    "custom_text_gene_list",
    "duo",
    "extraction",
    "genomic_intervals_collection",
    "max_variant",
    "pathology_test_version",
    "patient",
    "pedigree",
    "quad",
    "sample",
    "sample_gene_list",
    "specimen",
    "tissue_sample",
    "trio",
    "variant_collection",
]

# Node config keyed by the exporting server's own cohorts, samples and VCFs - the node arrives without those,
# so the rows are left out, and named in the export so the import can say what it didn't bring
NON_PORTABLE_NODE_CONFIG = {
    "cohortnodezygosityfilterscollection_set": "cohort zygosity filters",
    "nodevcffilter_set": "VCF filters",
    "samplenodesamplefilter_set": "per sample filters",
}

# Node PKs within the analysis being exported - remapped onto the new nodes after they're created
INTERNAL_NODE_REFERENCE_FIELDS = ["left_parent", "right_parent"]

# Runtime state - reload_analysis_nodes recalculates it at the end of the import
NODE_EXCLUDE = ["created", "modified", "version", "appearance_version", "ready", "valid", "count",
                "errors", "status", "load_seconds"]

ANALYSIS_EXCLUDE = ["id", "created", "modified", "user", "default_sort_by_column",
                    "custom_columns_collection", "tag_config_collection", "canonical_transcript_collection",
                    "visible", "annotation_version"]

# Exports written before the ZygosityRequirement choice (#2659) carry the boolean it replaced:
# {model: {old boolean field: new choice field}}
LEGACY_ZYGOSITY_REQUIREMENT_FIELDS = {
    "analysis.DuoNode": {"require_zygosity": "relative_zygosity_requirement"},
    "analysis.PedigreeNode": {"require_zygosity": "zygosity_requirement"},
    "analysis.QuadNode": {"require_parent_zygosity": "parent_zygosity_requirement",
                          "require_sibling_zygosity": "sibling_zygosity_requirement"},
    "analysis.TrioNode": {"require_zygosity": "parent_zygosity_requirement"},
}


def upgrade_legacy_node_fields(model_name: str, data: dict):
    for old_field, new_field in LEGACY_ZYGOSITY_REQUIREMENT_FIELDS.get(model_name, {}).items():
        if old_field in data:
            required = data.pop(old_field)
            data[new_field] = ZygosityRequirement.CALL_REQUIRED if required else ZygosityRequirement.ALLOW_NO_CALL


def analysis_export_to_dict(analysis: Analysis) -> dict:
    node_serializers = AnalysisNodeSerializer.get_node_serializers()

    analysis_serializer = AnalysisSerializer(analysis, exclude=ANALYSIS_EXCLUDE)

    analysis_data = {
        "analysis": {"model": analysis._meta.label, "pk": analysis.pk, "fields": dict(analysis_serializer.data)}
    }

    nodes = []
    for node in analysis.analysisnode_set.all().select_subclasses():
        model_name = node._meta.label
        serializer = node_serializers.get(model_name, AnalysisNodeSerializer)
        fields = dict(serializer(node, exclude=NODE_EXCLUDE + NON_PORTABLE_NODE_FIELDS).data)
        node_record = {"model": model_name, "pk": node.pk, "fields": dict(fields)}
        if not_exported := [description for accessor, description in NON_PORTABLE_NODE_CONFIG.items()
                            if hasattr(node, accessor) and getattr(node, accessor).exists()]:
            node_record["not_exported"] = not_exported
        nodes.append(node_record)

    analysis_data["nodes"] = nodes
    analysis_data["edges"] = json.loads(serialize('json', AnalysisEdge.objects.filter(child__analysis=analysis)))
    if node_count_config := AnalysisNodeCountConfiguration.objects.filter(analysis=analysis).first():
        node_counts_qs = node_count_config.analysisnodecountconfigrecord_set.order_by("sort_order")
        analysis_data["node_counts"] = list(node_counts_qs.values_list("node_count_type", flat=True))
    return analysis_data


def analysis_export_to_file(analysis: Analysis, file):
    analysis_data = analysis_export_to_dict(analysis)
    json.dump(analysis_data, file)


def analysis_import(user: User, genome_build: Optional[GenomeBuild], filename,
                    annotation_version: AnnotationVersion = None, replace: bool = True) -> Optional[Analysis]:
    with open(filename, encoding="utf-8") as f:
        analysis_json = json.loads(f.read())

    if not replace:
        name = analysis_json["analysis"]["fields"]["name"]
        if Analysis.objects.filter(name=name, template_type=AnalysisTemplateType.TEMPLATE).exists():
            return None

    import_warnings = []
    analysis = _create_analysis(user, genome_build, annotation_version, analysis_json, import_warnings)
    if import_warnings:
        details = f"Analysis {analysis.pk} '{analysis.name}' imported without:\n" + "\n".join(import_warnings)
        create_event(user, "analysis_import", details, filename=filename, severity=LogLevel.WARNING)
    reload_analysis_nodes(analysis.pk)
    return analysis


@transaction.atomic
def _create_analysis(user: User, genome_build: Optional[GenomeBuild],
                     annotation_version: Optional[AnnotationVersion], analysis_json: dict,
                     import_warnings: list[str]) -> Analysis:
    analysis_kwargs = analysis_json["analysis"]["fields"]
    exported_genome_build_name = analysis_kwargs.pop("genome_build", None)
    if genome_build is None:
        genome_build = GenomeBuild.get_name_or_alias(exported_genome_build_name)

    for field in ANALYSIS_EXCLUDE:
        analysis_kwargs.pop(field, None)  # Files written before these were excluded still carry them
    analysis_kwargs["template_type"] = None  # An import is a normal analysis - the user can make it a template

    analysis = Analysis(genome_build=genome_build, **analysis_kwargs)
    analysis.set_defaults_and_save(user)
    if annotation_version:
        analysis.annotation_version = annotation_version
        analysis.save()
    if (node_counts := analysis_json.get("node_counts")) is not None:
        analysis.set_node_count_types(node_counts)  # The file's configuration beats the user's defaults
    assign_permission_to_user_and_groups(user, analysis)

    node_serializers = AnalysisNodeSerializer.get_node_serializers()
    serializer_context = {'user': user, "import_warnings": import_warnings}

    logging.info("Creating analysis: %s", analysis.pk)
    old_new_map = {}
    old_internal_references = {}

    for node_record in analysis_json["nodes"]:
        model_name = node_record["model"]
        data = node_record["fields"]
        old_pk = data.pop("id")
        data["analysis"] = analysis.pk

        for field in NON_PORTABLE_NODE_FIELDS:
            data.pop(field, None)
        upgrade_legacy_node_fields(model_name, data)
        internal_references = {}
        for field in INTERNAL_NODE_REFERENCE_FIELDS:
            if old_reference_pk := data.pop(field, None):
                internal_references[field] = old_reference_pk
        if internal_references:
            old_internal_references[old_pk] = internal_references

        serializer = node_serializers[model_name]
        s = serializer(data=data, context=serializer_context)
        if not s.is_valid():
            raise ValueError(f"{model_name}: {s.errors}")
        node = s.save()
        if not_exported := node_record.get("not_exported"):
            # The exported name - a source node's is generated from the sample etc. it no longer has
            import_warnings.append(f"{node._meta.verbose_name} '{data.get('name', '')}': {', '.join(not_exported)} "
                                   "(they refer to the exporting server's own data)")
        old_new_map[old_pk] = node.pk

    for edge_record in analysis_json["edges"]:
        edge_fields = edge_record["fields"]
        old_parent_id = edge_fields["parent"]
        old_child_id = edge_fields["child"]

        parent_id = old_new_map.get(old_parent_id)
        child_id = old_new_map.get(old_child_id)
        if not (parent_id and child_id):
            raise ValueError(f"Edge parent: {old_parent_id} to child: {old_child_id} - node(s) missing from file")
        AnalysisEdge.objects.create(parent_id=parent_id, child_id=child_id)

    for old_pk, internal_references in old_internal_references.items():
        new_references = {}
        for field, old_reference_pk in internal_references.items():
            if new_reference_pk := old_new_map.get(old_reference_pk):
                new_references[f"{field}_id"] = new_reference_pk
            else:
                raise ValueError(f"Node {old_pk} {field}: {old_reference_pk} - node missing from file")
        VennNode.objects.filter(pk=old_new_map[old_pk]).update(**new_references)

    return analysis
