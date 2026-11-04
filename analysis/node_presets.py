"""
ACMG criterion presets for the Effect, Population and ClinVar node editors (#2605).

A preset is a named criterion (PVS1, PP3, PM2 ...) whose values fill the editor's form client side for the user to
adjust and save. Nothing is stored on the node, so its name, chips and method summary describe what it really
filters on. The Effect and ClinVar cutoffs are published calibrations, the ones the Pathogenicity Thresholds page
lists (`annotation/pathogenicity_predictions.py`), so they are set in one place. The Population (PM2) lines are lab
policy that varies with mode of inheritance, penetrance and severity, so they are `PopulationNodePreset` rows, seeded
by migration and edited in Django admin. Presets look for pathogenic evidence only (keep candidates, drop benign
evidence), which is what an analysis wants almost every time; the criteria on one node pull in different directions
(PM2 vs BA1, PVS1 vs PP3), so each preset is one criterion.

Entry points: `damage_node_presets`, `population_node_presets` and `clinvar_node_presets`, each taking the node being
edited (what it offers depends on the annotation version and the data imported) and returning the NodePresets its
editor lists. `analysis/forms/forms_nodes.py:BaseNodeForm.get_preset_data` turns each into a complete form state for
`node_presets.js`.
"""
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from analysis.enums import GroupOperation
from analysis.models.nodes.filters.population_node import PopulationNode, PopulationNodePreset
from annotation.enums import PathogenicityImpact
from annotation.pathogenicity_predictions import (
    GENE_LOEUF_CUTOFF,
    GENE_MIS_Z_CUTOFF,
    PathogenicityTool,
    RawScoreDirection,
)
from snpdb.enums import AlleleOriginFilterDefault

if TYPE_CHECKING:
    from analysis.models.nodes.filters.clinvar_node import ClinVarNode
    from analysis.models.nodes.filters.damage_node import DamageNode


@dataclass(frozen=True)
class NodePreset:
    key: str
    label: str  # the dropdown option, leading with the criterion
    description: str  # what the preset sets and why, shown once it has been applied
    values: dict[str, Any]  # form field -> value; every other field goes back to the node's default


def _pp3_tools(node: "DamageNode") -> list[PathogenicityTool]:
    """ Tools this annotation version scores whose cutoff is ClinGen's PP3-supporting line. A tool that runs the
        other way (popEVE) or is signed (PromoterAI) has no such cutoff, so a higher=worse check is enough """
    return [t for t in node.get_raw_score_tools()
            if t.clingen_calibrated and t.raw_pathogenic_threshold is not None
            and t.raw_direction == RawScoreDirection.HIGHER]


def damage_node_presets(node: "DamageNode") -> list[NodePreset]:
    presets = []

    pvs1 = {"impact_min": PathogenicityImpact.HIGH, "impact_required": True}
    pvs1_description = ("Impact HIGH (stop gained, frameshift, canonical splice site, start lost), required. ")
    if node.has_gene_constraint:
        pvs1 |= {"gene_loeuf_max": GENE_LOEUF_CUTOFF, "gene_loeuf_required": True, "gene_loeuf_allow_null": True}
        pvs1_description += (f"gnomAD LOEUF ≤ {GENE_LOEUF_CUTOFF:g}, required - loss of function is a known "
                             f"mechanism for the gene - keeping genes with no constraint score. ")
    pvs1_description += ("The PTC NMD filter is left off: it only applies to frameshifts and stop gains whose new "
                         "stop was located, so set it to narrow to those.")
    presets.append(NodePreset("pvs1", "PVS1 - null variant candidates", pvs1_description, pvs1))

    if pp3_tools := _pp3_tools(node):
        pp3 = {"impact_min": PathogenicityImpact.MODERATE, "impact_required": True}
        for tool in pp3_tools:
            pp3[tool.node_threshold_field] = tool.raw_pathogenic_threshold
            pp3[f"{tool.raw_field}_allow_null"] = True
        cutoffs = ", ".join(f"{t.name} ≥ {t.raw_pathogenic_threshold:g}" for t in pp3_tools)
        presets.append(NodePreset(
            "pp3", "PP3 - missense with computational support",
            f"Impact MODERATE or above, required, then any ClinGen-calibrated predictor at its PP3-supporting "
            f"cutoff: {cutoffs}. A variant a tool has no score for is kept when another tool supports it.",
            pp3))

    if node.has_gene_constraint:
        presets.append(NodePreset(
            "pp2", "PP2 - missense in a constrained gene",
            f"Impact MODERATE or above, required (missense and in-frame indels), with gnomAD missense Z ≥ "
            f"{GENE_MIS_Z_CUTOFF:g}, required - the gene has fewer missense variants than expected, so missense "
            f"is a common mechanism of disease. Genes with no constraint score are dropped.",
            {"impact_min": PathogenicityImpact.MODERATE, "impact_required": True,
             "gene_mis_z_min": GENE_MIS_Z_CUTOFF, "gene_mis_z_required": True, "gene_mis_z_allow_null": False}))
    return presets


def _population_preset_summary(values: dict[str, Any]) -> str:
    databases = [label for field, label in PopulationNode.POPULATION_DATABASE_SHORT_LABELS.items() if values.get(field)]
    parts = []
    if databases and values["percent"] != PopulationNode.EVERYTHING:
        parts.append(f"{' and '.join(databases)} ≤ {values['percent']:g}%")
    for field, label in PopulationNode.GNOMAD_COUNT_LABELS.items():
        if (max_count := values.get(field)) is not None:
            parts.append(f"{label} ≤ {max_count}")
    summary = "; ".join(parts) or "no population filter"
    if values["keep_internally_classified_pathogenic"]:
        summary += ", keeping variants this database has classified (likely) pathogenic"
    return f"Sets {summary}."


def population_node_presets(node: PopulationNode) -> list[NodePreset]:
    presets = []
    for preset in PopulationNodePreset.objects.filter(enabled=True):
        values = {field: getattr(preset, field) for field in PopulationNodePreset.NODE_FIELDS}
        values["group_operation"] = GroupOperation.ANY
        if not node.has_filtering_allele_frequency:  # the editor doesn't draw them, so they must not filter
            values["gnomad_fafmax_faf95_max"] = values["gnomad_fafmax_faf99_max"] = False
        description = " ".join(filter(None, [preset.description, _population_preset_summary(values)]))
        presets.append(NodePreset(f"population_{preset.pk}", preset.name, description, values))
    return presets


def clinvar_node_presets(_node: "ClinVarNode") -> list[NodePreset]:
    return [NodePreset(
        "pp5", "PP5 - reputable source reports pathogenic",
        "Germline Pathogenic or Likely pathogenic with at least 2 review stars (multiple submitters with no "
        "conflicts, an expert panel or a practice guideline). Records with conflicting interpretations are not "
        "matched.",
        {"allele_origin": AlleleOriginFilterDefault.GERMLINE,
         "germline_pathogenic": True, "germline_likely_pathogenic": True, "germline_uncertain": False,
         "germline_likely_benign": False, "germline_benign": False, "germline_other": False,
         "stars_min": 2, "conflicting": False})]
