"""
`vg inspect analysis <pk>`: an Analysis and its node DAG - build, annotation version, template
type, lock, then every node in topological order with class, name, version, status, count, the
last error and whether a NodeVersion / cache exists, plus edges, tags and template provenance.
"""
from typing import Any

from analysis.models.models_analysis import Analysis, AnalysisTemplateRun
from analysis.models.models_variant_tag import VariantTag
from analysis.models.nodes.analysis_node import (
    AnalysisEdge,
    AnalysisNode,
    NodeCache,
    NodeTask,
    NodeVersion,
)
from analysis.models.nodes.node_utils import get_toposorted_nodes
from library.vg.inspect import capped
from library.vg.inspect.common import user_summary

NODE_CAP = 60


def load(key: str) -> Analysis:
    try:
        return Analysis.objects.select_related("genome_build", "user", "annotation_version").get(pk=int(key))
    except (ValueError, Analysis.DoesNotExist) as e:
        raise LookupError(f"No Analysis with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    analysis = load(key)
    data: dict[str, Any] = {
        "id": analysis.pk,
        "name": analysis.name,
        "user": user_summary(analysis.user) if analysis.user_id else None,
        "build": analysis.genome_build_id,
        "annotation_version": analysis.annotation_version_id,
        "template_type": analysis.get_template_type_display() if analysis.template_type else None,
        "version": analysis.version,
        "visible": analysis.visible,
        "locked": analysis.is_locked,
        "created": analysis.created.date().isoformat(),
        "modified": analysis.modified.date().isoformat(),
        "errors": analysis.get_errors(),
        "warnings": analysis.get_warnings(),
        "url": analysis.get_absolute_url(),
    }
    run = AnalysisTemplateRun.objects.filter(analysis=analysis).select_related("template_version__template").first()
    if run:
        data["from_template"] = {"template": run.template_version.template.name, "version": run.template_version.version}
    nodes_qs = AnalysisNode.objects.filter(analysis=analysis).select_subclasses()
    nodes = list(nodes_qs)
    data["node_count"] = len(nodes)
    if depth >= 2:
        ordered = [node for group in get_toposorted_nodes(nodes_qs) for node in group]
        data["nodes"] = capped(ordered, _node, cap=NODE_CAP)
        data["edges"] = capped(AnalysisEdge.objects.filter(parent__analysis=analysis).order_by("parent_id", "child_id"),
                               lambda e: f"{e.parent_id} -> {e.child_id}", cap=NODE_CAP * 2)
        data["tags"] = capped(VariantTag.objects.filter(analysis=analysis).select_related("tag").order_by("-pk"),
                              lambda t: {"tag": t.tag_id, "variant": t.variant_id, "node": t.node_id})
    return data


def _node(node: AnalysisNode) -> dict[str, Any]:
    info: dict[str, Any] = {"id": node.pk, "class": node.get_class_name(), "name": node.name or node.get_name_or_identifier(),
                            "version": node.version, "status": node.get_status_display(), "count": node.count}
    if node.errors:
        info["error"] = node.errors.strip().splitlines()[-1][:160]
    node_version = NodeVersion.objects.filter(node=node, version=node.version).first()
    if node_version:
        info["counts"] = node_version.load_data.get("counts") if node_version.load_data else None
        info["cache"] = NodeCache.objects.filter(node_version=node_version).exists()
        task = NodeTask.objects.filter(node_version=node_version).first()
        if task:
            info["task"] = task.celery_task or "(lease)"
    else:
        info["node_version"] = "missing"
    return info
