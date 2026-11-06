import itertools
import operator
from functools import cached_property, reduce
from typing import NamedTuple, Optional

from auditlog.registry import auditlog
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models.query_utils import Q

from analysis.enums import (
    ALoFTPredictionOptions,
    PTCNMDFilter,
    SpliceVaultOutOfFrame,
    UTRHighImpact,
)
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_display import NodeChip, NodeIcon
from annotation.annotsv_columns import ANNOTSV_ARTEFACT_REGION_COLUMNS
from annotation.effect_filters import (
    gene_constraint_q,
    splicevault_out_of_frame_q,
    utr_high_impact_q,
)
from annotation.enums import (
    ALoFTPrediction,
    AlphaMissensePrediction,
    ClinPredPrediction,
    EVEClass,
    MetaRNNPrediction,
    Pathogenicity,
    PathogenicityImpact,
    PrimateAIPrediction,
    VariantAnnotationPipelineType,
)
from annotation.models.models import VariantAnnotation
from annotation.models.models_annotation_run import AnnotationPipelineVersion
from annotation.models.models_annotation_version import AnnotationVersion
from annotation.models.models_gene_annotation import GeneAnnotation
from annotation.pathogenicity_predictions import (
    TOOLS,
    TOOLS_BY_PRED_FIELD,
    PathogenicityTool,
    RawScoreDirection,
)
from genes.models.models_gnomad_gene_constraint import GnomADGeneConstraint
from snpdb.enums import BuiltInFilters
from snpdb.variant_filters import (
    get_variant_type_labels,
    get_variant_types_q,
    restricts_variant_types,
)


class _ColumnFilter(NamedTuple):
    """ A DamageNode filter over one annotation column """
    column: str  # VariantAnnotation column the filter reads
    value_field: str
    required_field: str
    prediction: bool = True  # a per-tool pathogenicity prediction, rather than a LoF flag
    rankscore: bool = False
    numeric: bool = False  # a threshold, where 0 is a value the user chose (the editor's enable box clears to None)

    def is_set(self, node) -> bool:
        value = getattr(node, self.value_field)
        return value is not None if self.numeric else bool(value)


_COLUMN_FILTERS = (
    *(_ColumnFilter(t.rankscore_field, f"{t.rankscore_field}_min", f"{t.rankscore_field}_required",
                    rankscore=True, numeric=True)
      for t in TOOLS if t.rankscore_field),
    *(_ColumnFilter(t.raw_field, t.node_threshold_field, f"{t.raw_field}_required", numeric=True)
      for t in TOOLS if t.raw_field),
    *(_ColumnFilter(t.pred_field, t.pred_field, f"{t.pred_field}_required") for t in TOOLS if t.pred_field),
    _ColumnFilter("nmd_escaping_variant", "nmd_escaping_variant", "nmd_escaping_variant_required", prediction=False),
    _ColumnFilter("aloft_pred", "aloft", "aloft_required", prediction=False),
    # AlphaGenome splicing is an OR term inside the splicing filter, so it shares splice_required
    _ColumnFilter("alphagenome_splicing", "alphagenome_splicing_min", "splice_required", prediction=False),
    _ColumnFilter("splicevault_out_of_frame_events", "splicevault_out_of_frame", "splicevault_out_of_frame_required",
                  prediction=False),
    _ColumnFilter("utr_annotator_annotation", "utr_high_impact", "utr_high_impact_required", prediction=False),
)

# Added to annotation versions already in use (#2426, #2433): null on variants annotated before the switch-on until
# VariantAnnotationVersion.backfilled_utr_avi_splicevault, so filtering on them earlier would drop every such variant
_PLUGIN_BACKFILL_COLUMNS = frozenset({"avi_phred", "alphagenome_splicing", "splicevault_out_of_frame_events",
                                      "utr_annotator_annotation"})

# Composite grid columns (VariantGridColumn) showing what the AnnotSV filters read
ANNOTSV_ACMG_GRID_COLUMN = "annotsv_acmg"
ANNOTSV_REGION_GRID_COLUMN = "annotsv_region"


def get_omim_q(annotation_version: AnnotationVersion) -> Q:
    """ Genes with OMIM terms, as a subquery against the small GeneAnnotation table - going through
        'variantannotation__gene__geneannotation' joins the huge VariantAnnotation table through Gene """
    gene_annotation_qs = GeneAnnotation.objects.filter(omim_terms__isnull=False)
    if gene_annotation_version_id := annotation_version.gene_annotation_version_id:
        gene_annotation_qs = gene_annotation_qs.filter(version_id=gene_annotation_version_id)
    return Q(variantannotation__gene__in=gene_annotation_qs.values_list("gene_id", flat=True))


class DamageNode(AnalysisNode):
    """ This is called 'EffectNode' in analysis """
    COUNT_LABEL_FIELDS = {
        BuiltInFilters.IMPACT_HIGH_OR_MODERATE: {"impact_min": PathogenicityImpact.MODERATE},
        BuiltInFilters.COSMIC: {"in_cosmic": True},
        BuiltInFilters.OMIM: {"omim_phenotype": True},
    }

    # The All Variants grid's variant types (snpdb.variant_filters), read off the stored alt rather than VEP's
    # variant_class, which calls a 1 Mb <DEL> and a 1 bp deletion the same thing (#2674). A type restriction
    # restricts - it goes in and_filters rather than the scoring OR pool, so it has no _required flag.
    # Empty, or every type, means no restriction.
    variant_types = ArrayField(models.CharField(max_length=20), default=list, blank=True)

    impact_min = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    impact_required = models.BooleanField(default=False)

    splice_min = models.FloatField(null=True, blank=True)
    splice_required = models.BooleanField(default=False)
    splice_allow_null = models.BooleanField(default=True)
    # GRCh38 SNVs, licence-restricted data (#2561). Its own threshold (0 to ~3) but part of the splicing filter
    alphagenome_splicing_min = models.FloatField(null=True, blank=True)

    cosmic_count_min = models.IntegerField(null=True, blank=True)
    cosmic_count_required = models.BooleanField(default=False)

    # "Has a COSMIC id", separate from cosmic_count_min: the id comes from VEP's Existing_variation on every
    # deployment, the count only from the licensed COSMIC custom VCF
    in_cosmic = models.BooleanField(default=False)
    in_cosmic_required = models.BooleanField(default=False)
    # The variant's gene has OMIM phenotype terms (GeneAnnotation.omim_terms for the analysis's gene annotation version)
    omim_phenotype = models.BooleanField(default=False)
    omim_phenotype_required = models.BooleanField(default=False)

    protein_domain = models.BooleanField(default=False)
    protein_domain_required = models.BooleanField(default=False)

    published = models.BooleanField(default=False)
    published_required = models.BooleanField(default=False)

    damage_predictions_min = models.IntegerField(null=True, blank=True)
    damage_predictions_required = models.BooleanField(default=False)
    damage_predictions_allow_null = models.BooleanField(default=True)

    # Columns v1
    cadd_score_min = models.IntegerField(null=True, blank=True)
    cadd_score_required = models.BooleanField(default=False)
    cadd_score_allow_null = models.BooleanField(default=True)

    revel_score_min = models.FloatField(null=True, blank=True)
    revel_score_required = models.BooleanField(default=False)
    revel_score_allow_null = models.BooleanField(default=True)

    # Rankscores added in columns v2, alpha missense in v3
    alphamissense_rankscore_min = models.FloatField(null=True, blank=True)
    alphamissense_rankscore_required = models.BooleanField(default=False)
    alphamissense_rankscore_allow_null = models.BooleanField(default=True)

    bayesdel_noaf_rankscore_min = models.FloatField(null=True, blank=True)
    bayesdel_noaf_rankscore_required = models.BooleanField(default=False)
    bayesdel_noaf_rankscore_allow_null = models.BooleanField(default=True)

    cadd_raw_rankscore_min = models.FloatField(null=True, blank=True)
    cadd_raw_rankscore_required = models.BooleanField(default=False)
    cadd_raw_rankscore_allow_null = models.BooleanField(default=True)

    clinpred_rankscore_min = models.FloatField(null=True, blank=True)
    clinpred_rankscore_required = models.BooleanField(default=False)
    clinpred_rankscore_allow_null = models.BooleanField(default=True)

    metalr_rankscore_min = models.FloatField(null=True, blank=True)
    metalr_rankscore_required = models.BooleanField(default=False)
    metalr_rankscore_allow_null = models.BooleanField(default=True)

    revel_rankscore_min = models.FloatField(null=True, blank=True)
    revel_rankscore_required = models.BooleanField(default=False)
    revel_rankscore_allow_null = models.BooleanField(default=True)

    vest4_rankscore_min = models.FloatField(null=True, blank=True)
    vest4_rankscore_required = models.BooleanField(default=False)
    vest4_rankscore_allow_null = models.BooleanField(default=True)

    # Raw scores added in columns v4 (#1543)
    alphamissense_score_min = models.FloatField(null=True, blank=True)
    alphamissense_score_required = models.BooleanField(default=False)
    alphamissense_score_allow_null = models.BooleanField(default=True)

    bayesdel_noaf_score_min = models.FloatField(null=True, blank=True)
    bayesdel_noaf_score_required = models.BooleanField(default=False)
    bayesdel_noaf_score_allow_null = models.BooleanField(default=True)

    # cadd_phred_min reuses the user-facing CADD scale; cadd_raw is not exposed.
    cadd_phred_min = models.FloatField(null=True, blank=True)
    cadd_phred_required = models.BooleanField(default=False)
    cadd_phred_allow_null = models.BooleanField(default=True)

    clinpred_score_min = models.FloatField(null=True, blank=True)
    clinpred_score_required = models.BooleanField(default=False)
    clinpred_score_allow_null = models.BooleanField(default=True)

    # The VEP 116 plugin scores below are GRCh38 + columns_version >= 5, so unlike the dbNSFP scores
    # their sliders only show on annotation versions that populate them (#1808)
    eve_score_min = models.FloatField(null=True, blank=True)
    eve_score_required = models.BooleanField(default=False)
    eve_score_allow_null = models.BooleanField(default=True)

    # popEVE is a log-likelihood ratio - more negative is more damaging, hence _max not _min
    popeve_score_max = models.FloatField(null=True, blank=True)
    popeve_score_required = models.BooleanField(default=False)
    popeve_score_allow_null = models.BooleanField(default=True)

    # ProtVar ddG (kcal/mol) and PromoterAI aren't missense pathogenicity predictions - they filter
    # only, and stay out of predictions_num_pathogenic. PromoterAI's score is signed, so it matches
    # on |score|.
    protvar_stability_min = models.FloatField(null=True, blank=True)
    protvar_stability_required = models.BooleanField(default=False)
    protvar_stability_allow_null = models.BooleanField(default=True)

    promoter_ai_score_min = models.FloatField(null=True, blank=True)
    promoter_ai_score_required = models.BooleanField(default=False)
    promoter_ai_score_allow_null = models.BooleanField(default=True)

    # AVI PHRED is a genome-wide rank with no published cutoff (#2561) - filter only, like PromoterAI
    avi_phred_min = models.FloatField(null=True, blank=True)
    avi_phred_required = models.BooleanField(default=False)
    avi_phred_allow_null = models.BooleanField(default=True)

    metarnn_score_min = models.FloatField(null=True, blank=True)
    metarnn_score_required = models.BooleanField(default=False)
    metarnn_score_allow_null = models.BooleanField(default=True)

    mpc_score_min = models.FloatField(null=True, blank=True)
    mpc_score_required = models.BooleanField(default=False)
    mpc_score_allow_null = models.BooleanField(default=True)

    mutpred2_score_min = models.FloatField(null=True, blank=True)
    mutpred2_score_required = models.BooleanField(default=False)
    mutpred2_score_allow_null = models.BooleanField(default=True)

    primateai_score_min = models.FloatField(null=True, blank=True)
    primateai_score_required = models.BooleanField(default=False)
    primateai_score_allow_null = models.BooleanField(default=True)

    varity_er_score_min = models.FloatField(null=True, blank=True)
    varity_er_score_required = models.BooleanField(default=False)
    varity_er_score_allow_null = models.BooleanField(default=True)

    varity_r_score_min = models.FloatField(null=True, blank=True)
    varity_r_score_required = models.BooleanField(default=False)
    varity_r_score_allow_null = models.BooleanField(default=True)

    vest4_score_min = models.FloatField(null=True, blank=True)
    vest4_score_required = models.BooleanField(default=False)
    vest4_score_allow_null = models.BooleanField(default=True)

    alphamissense_pred = models.CharField(max_length=1, choices=AlphaMissensePrediction.choices, null=True, blank=True)
    alphamissense_pred_required = models.BooleanField(default=False)
    alphamissense_pred_allow_null = models.BooleanField(default=True)

    clinpred_pred = models.CharField(max_length=1, choices=ClinPredPrediction.CHOICES, null=True, blank=True)
    clinpred_pred_required = models.BooleanField(default=False)
    clinpred_pred_allow_null = models.BooleanField(default=True)

    metarnn_pred = models.CharField(max_length=1, choices=MetaRNNPrediction.CHOICES, null=True, blank=True)
    metarnn_pred_required = models.BooleanField(default=False)
    metarnn_pred_allow_null = models.BooleanField(default=True)

    primateai_pred = models.CharField(max_length=1, choices=PrimateAIPrediction.CHOICES, null=True, blank=True)
    primateai_pred_required = models.BooleanField(default=False)
    primateai_pred_allow_null = models.BooleanField(default=True)

    # EVE's class is stored as a word, not the dbNSFP single-char code
    eve_class = models.CharField(max_length=10, choices=EVEClass.choices, null=True, blank=True)
    eve_class_required = models.BooleanField(default=False)
    eve_class_allow_null = models.BooleanField(default=True)

    nmd_escaping_variant = models.BooleanField(default=False)
    nmd_escaping_variant_required = models.BooleanField(default=False)

    # Filters on nmd_escape_status (columns v5, #579, #2596)
    ptc_nmd = models.CharField(max_length=1, choices=PTCNMDFilter.choices, null=True, blank=True)
    ptc_nmd_required = models.BooleanField(default=False)

    aloft = models.CharField(max_length=1, choices=ALoFTPredictionOptions.choices, null=True, blank=True)
    aloft_required = models.BooleanField(default=False)
    aloft_allow_null = models.BooleanField(default=True)

    # #2561 - categorical plugin annotations, see annotation/effect_filters.py
    splicevault_out_of_frame = models.CharField(max_length=1, choices=SpliceVaultOutOfFrame.choices,
                                                null=True, blank=True)
    splicevault_out_of_frame_required = models.BooleanField(default=False)
    splicevault_out_of_frame_allow_null = models.BooleanField(default=True)

    utr_high_impact = models.CharField(max_length=1, choices=UTRHighImpact.choices, null=True, blank=True)
    utr_high_impact_required = models.BooleanField(default=False)
    utr_high_impact_allow_null = models.BooleanField(default=True)

    # #2600 - gnomAD gene constraint of the variant's gene (genes.GnomADGeneConstraint, not a VEP column).
    # LOEUF: PVS1's "LoF is a known mechanism"; missense Z: PP2 / BP1
    gene_loeuf_max = models.FloatField(null=True, blank=True)
    gene_loeuf_required = models.BooleanField(default=False)
    gene_loeuf_allow_null = models.BooleanField(default=True)

    gene_mis_z_min = models.FloatField(null=True, blank=True)
    gene_mis_z_required = models.BooleanField(default=False)
    gene_mis_z_allow_null = models.BooleanField(default=True)

    # PM4's "non-repeat region" - a restriction like variant_types, so it has no _required flag
    outside_repeat_region = models.BooleanField(default=False)

    # #2701 - AnnotSV's SV annotations, null on every small variant. Required + allow null reads "SVs must be LP/P,
    # small variants untouched"; in the OR pool, "damaging by any evidence, AnnotSV included"
    annotsv_acmg_class_min = models.IntegerField(choices=Pathogenicity.choices, null=True, blank=True)
    annotsv_acmg_class_required = models.BooleanField(default=False)
    annotsv_acmg_class_allow_null = models.BooleanField(default=True)

    annotsv_acmg_score_min = models.FloatField(null=True, blank=True)
    annotsv_acmg_score_required = models.BooleanField(default=False)
    annotsv_acmg_score_allow_null = models.BooleanField(default=True)

    # No SV breakpoint in ANNOTSV_ARTEFACT_REGION_COLUMNS - a restriction like outside_repeat_region
    exclude_sv_artefact_regions = models.BooleanField(default=False)

    def get_raw_score_tools(self) -> list[PathogenicityTool]:
        """ Tools whose raw-score slider this analysis's annotation version actually populates -
            gated off visible_columns so the editor can't offer a filter over an empty column. """
        visible_columns = self.analysis.annotation_version.variant_annotation_version.visible_columns
        return [t for t in TOOLS if t.raw_field in visible_columns]

    def get_pred_tools(self) -> list[PathogenicityTool]:
        """ Tools offering a categorical prediction dropdown on this annotation version """
        visible_columns = self.analysis.annotation_version.variant_annotation_version.visible_columns
        return [t for t in TOOLS if t.pred_field in visible_columns]

    def _column_filters(self) -> list[_ColumnFilter]:
        """ The per-column filters this analysis's annotation version offers - gated off visible_columns
            (VEP_COLUMNS), so the version and build decide, as they do for the editor's sliders. """
        if self.columns_version == 1:
            # cadd_phred / revel_score are populated again from v4, where the node filters them as raw scores
            return [_ColumnFilter("cadd_phred", "cadd_score_min", "cadd_score_required", numeric=True),
                    _ColumnFilter("revel_score", "revel_score_min", "revel_score_required", numeric=True)]
        visible_columns = self.analysis.annotation_version.variant_annotation_version.visible_columns
        column_filters = [cf for cf in _COLUMN_FILTERS if cf.column in visible_columns]
        # nmd_escape_status is calculated by VariantGrid, not a VEP column
        if self.columns_version >= 5:
            column_filters.append(_ColumnFilter("nmd_escape_status", "ptc_nmd", "ptc_nmd_required",
                                                prediction=False))
        return column_filters

    def has_alphagenome_splicing(self) -> bool:
        return "alphagenome_splicing" in self.analysis.annotation_version.variant_annotation_version.visible_columns

    def has_splicevault(self) -> bool:
        return ("splicevault_out_of_frame_events"
                in self.analysis.annotation_version.variant_annotation_version.visible_columns)

    def has_utr_annotator(self) -> bool:
        return "utr_annotator_annotation" in self.analysis.annotation_version.variant_annotation_version.visible_columns

    @cached_property
    def has_gene_constraint(self) -> bool:
        """ The gnomAD gene constraint is a CachedWebResource import, so a deployment may not have it """
        return GnomADGeneConstraint.objects.exists()

    def has_repeat_masker(self) -> bool:
        return "repeat_masker" in self.analysis.annotation_version.variant_annotation_version.visible_columns

    def _gene_constraint_filters(self) -> list[tuple[str, str, Optional[float], bool, bool]]:
        """ (constraint field, lookup, threshold, required, allow_null) for each gene constraint filter set """
        if not self.has_gene_constraint:
            return []
        filters = [("lof_oe_ci_upper", "lte", self.gene_loeuf_max, self.gene_loeuf_required,
                    self.gene_loeuf_allow_null),
                   ("mis_z_score", "gte", self.gene_mis_z_min, self.gene_mis_z_required, self.gene_mis_z_allow_null)]
        return [f for f in filters if f[2] is not None]

    def _applies_outside_repeat_region(self) -> bool:
        return self.outside_repeat_region and self.has_repeat_masker()

    @cached_property
    def has_annotsv(self) -> bool:
        """ AnnotSV is its own pipeline rather than VEP columns, so visible_columns doesn't cover it """
        genome_build = self.analysis.genome_build
        return AnnotationPipelineVersion.get_active(VariantAnnotationPipelineType.ANNOTSV, genome_build) is not None

    def _annotsv_filters_set(self) -> bool:
        return (self.annotsv_acmg_class_min is not None or self.annotsv_acmg_score_min is not None
                or self.exclude_sv_artefact_regions)

    def _applies_annotsv(self) -> bool:
        """ Checks the node's own settings first, so a node without AnnotSV filters never looks up the pipeline """
        return self._annotsv_filters_set() and self.has_annotsv

    def _annotsv_threshold_filters(self) -> list[tuple[str, float, bool, bool]]:
        """ (column, threshold, required, allow_null) for each AnnotSV threshold set and applied """
        if not self._applies_annotsv():
            return []
        filters = [("annotsv_acmg_class", self.annotsv_acmg_class_min, self.annotsv_acmg_class_required,
                    self.annotsv_acmg_class_allow_null),
                   ("annotsv_acmg_score", self.annotsv_acmg_score_min, self.annotsv_acmg_score_required,
                    self.annotsv_acmg_score_allow_null)]
        return [f for f in filters if f[1] is not None]

    def _applies_exclude_sv_artefact_regions(self) -> bool:
        return self.exclude_sv_artefact_regions and self._applies_annotsv()

    def _filter_column_applied(self, column: str) -> bool:
        """ False where the column is licence-restricted data the version withholds (the editor hides its slider, so a
            saved threshold can't be seen or changed - #2571) or one the version is still waiting on a backfill for.
            The filter is skipped (and get_warnings says so) rather than filtering on data the user can't see, or
            treating a not-yet-imported null as "not damaging" """
        vav = self.analysis.annotation_version.variant_annotation_version
        if column in vav.withheld_licence_restricted_columns:
            return False
        return column not in _PLUGIN_BACKFILL_COLUMNS or vav.backfilled_utr_avi_splicevault

    def modifies_parents(self):
        all_versions = [self.impact_min, self.splice_min, self.cosmic_count_min, self.damage_predictions_min,
                        self.in_cosmic, self.omim_phenotype, self.protein_domain, self.published,
                        restricts_variant_types(self.variant_types), bool(self._gene_constraint_filters()),
                        self._applies_outside_repeat_region(), self._applies_annotsv()]
        return any(all_versions) or any(cf.is_set(self) for cf in self._column_filters())

    def has_required(self) -> bool:
        all_versions = [self.impact_required, self.splice_required, self.cosmic_count_required,
                        self.in_cosmic_required, self.omim_phenotype_required,
                        self.damage_predictions_required, self.protein_domain_required, self.published_required,
                        *(required for _, _, _, required, _ in self._gene_constraint_filters()),
                        *(required for _, _, required, _ in self._annotsv_threshold_filters())]
        # A rankscore's required checkbox is only on screen while the rankscore sliders are
        show_rankscores = self.show_rankscore_sliders()
        return any(all_versions) or any(getattr(self, cf.required_field) for cf in self._column_filters()
                                        if show_rankscores or not cf.rankscore)

    def has_individual_pathogenic_predictions(self) -> bool:
        return any(cf.is_set(self) for cf in self._column_filters() if cf.prediction)

    def has_rankscore_predictions(self) -> bool:
        """ v4 hides rankscore sliders by default - used to auto-reveal them when one is still set. """
        return any(cf.is_set(self) for cf in _COLUMN_FILTERS if cf.rankscore)

    def show_rankscore_sliders(self) -> bool:
        """ v2/v3 only have rankscores so always show them. At v4 they're legacy: offered when the
        deployment allows legacy rankscores, or this node already has one set. From v5 they're
        legacy-only - just the nodes that already carry one (saved filters stay visible/editable,
        and the filter is applied regardless). """
        if self.columns_version < 4:
            return True
        if self.columns_version >= 5:
            return self.has_rankscore_predictions()
        return settings.ANNOTATION_SHOW_LEGACY_RANKSCORES or self.has_rankscore_predictions()

    def damage_predictions_description(self) -> str:
        return self.analysis.annotation_version.variant_annotation_version.damage_predictions_description

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        vav = self.analysis.annotation_version.variant_annotation_version
        if self.damage_predictions_min is not None and not vav.backfilled_damage_counts:
            warnings.append(
                "damage_predictions_min filter is NOT being applied on this VariantAnnotationVersion: "
                "predictions_num_pathogenic / predictions_num_benign aggregates have not been backfilled. "
                "Set the analysis to the latest annotation version to use the filter."
            )
        if self.ptc_nmd and self.columns_version >= 5 and not vav.backfilled_ptc:
            warnings.append(
                "PTC NMD filter is NOT being applied on this VariantAnnotationVersion: "
                "nmd_escape_status has not been backfilled. "
                "Set the analysis to the latest annotation version to use the filter."
            )
        set_filters = [(tool.name, tool.raw_field) for tool in TOOLS
                       if tool.raw_field and getattr(self, tool.node_threshold_field) is not None]
        if self.alphagenome_splicing_min is not None:
            set_filters.append(("AlphaGenome splicing", "alphagenome_splicing"))
        if self.splicevault_out_of_frame:
            set_filters.append(("SpliceVault out of frame", "splicevault_out_of_frame_events"))
        if self.utr_high_impact:
            set_filters.append(("High-impact 5' UTR", "utr_annotator_annotation"))
        for name, column in set_filters:
            if column in vav.withheld_licence_restricted_columns:
                warnings.append(
                    f"{name} filter is NOT being applied on this VariantAnnotationVersion: {column} is "
                    "licence-restricted data this version doesn't offer (not annotated with it, or its licence "
                    "rules out a use this deployment makes)."
                )
            elif column not in vav.visible_columns:
                warnings.append(
                    f"{name} filter has no data on this VariantAnnotationVersion: "
                    f"{column} is not populated on {vav.genome_build} / VEP {vav.vep} / "
                    f"columns_version {vav.columns_version}."
                )
            elif not self._filter_column_applied(column):
                warnings.append(
                    f"{name} filter is NOT being applied on this VariantAnnotationVersion: {column} has not "
                    "been backfilled for variants annotated before the plugin was switched on."
                )
        if self.utr_high_impact and not vav.backfilled_utr_effects:
            warnings.append(
                "High-impact 5' UTR filter may miss variants on this VariantAnnotationVersion: variants annotated "
                "by Ensembl's UTRAnnotator kept only one of their 5' UTR effects until the UTR effects backfill."
            )

        if (self.gene_loeuf_max is not None or self.gene_mis_z_min is not None) and not self.has_gene_constraint:
            warnings.append(
                "Gene constraint filter is NOT being applied: this server has not imported the gnomAD gene "
                "constraint data."
            )
        if self.outside_repeat_region and not self.has_repeat_masker():
            warnings.append(
                f"Outside repeat region filter is NOT being applied: repeat_masker is not populated on "
                f"{vav.genome_build} / VEP {vav.vep} / columns_version {vav.columns_version}."
            )
        if self._annotsv_filters_set() and not self.has_annotsv:
            warnings.append(
                f"AnnotSV filters are NOT being applied: there is no active AnnotSV annotation pipeline for "
                f"{self.analysis.genome_build}."
            )

        if self.splice_min is not None and vav.uses_raw_spliceai:
            warnings.append(
                "SpliceAI scores on this VariantAnnotationVersion are from the raw precomputed file "
                f"({vav.spliceai}). Illumina recommends the masked file for variant interpretation: "
                "the raw file also reports strengthening of annotated splice sites and weakening of "
                "unannotated ones, which are typically much less pathogenic. See "
                "https://github.com/Illumina/SpliceAI"
            )
        return warnings

    @staticmethod
    def _raw_score_q(tool: PathogenicityTool, threshold: float) -> Q:
        """ Keeps the damaging side of the score - which side that is depends on the tool """
        field = f"variantannotation__{tool.raw_field}"
        match tool.raw_direction:
            case RawScoreDirection.LOWER:
                return Q(**{f"{field}__lte": threshold})
            case RawScoreDirection.MAGNITUDE:
                return Q(**{f"{field}__gte": abs(threshold)}) | Q(**{f"{field}__lte": -abs(threshold)})
            case _:
                return Q(**{f"{field}__gte": threshold})

    def _get_node_q_hash(self) -> str:
        return str(self._get_node_q())

    def _get_node_q(self) -> Optional[Q]:
        or_filters = []
        and_filters = []

        if (q_variant_types := get_variant_types_q(self.variant_types)) is not None:
            and_filters.append(q_variant_types)

        if self.impact_min is not None:
            q_impact = PathogenicityImpact.get_q(self.impact_min)
            if self.impact_required:
                and_filters.append(q_impact)
            else:
                or_filters.append(q_impact)

        vav = self.analysis.annotation_version.variant_annotation_version
        splicing_q_list = []
        if self.splice_min is not None:
            # [consequence contains 'splice' OR not null splice region] AND [variant class not SNV]
            q_splice_indels = Q(variantannotation__consequence__contains='splice') | Q(variantannotation__splice_region__isnull=False)
            q_splice_indels &= Q(variantannotation__variant_class__ne="SN")
            splicing_q_list.extend([
                q_splice_indels,
                Q(variantannotation__dbscsnv_ada_score__gte=self.splice_min),
                Q(variantannotation__dbscsnv_rf_score__gte=self.splice_min),
            ])
            if self.splice_required and self.splice_allow_null:
                splicing_q_list.extend([
                    Q(variantannotation__dbscsnv_ada_score__isnull=True),
                    Q(variantannotation__dbscsnv_rf_score__isnull=True),
                ])

            # SpliceAI: use the optimised spliceai_max_ds column when the VAV has
            # been backfilled, otherwise fall back to the per-DS-field loop
            # (pre-c6f3c4e6f). See VariantAnnotationVersion.backfilled_spliceai_max_ds.
            if vav.backfilled_spliceai_max_ds:
                splicing_q_list.append(Q(variantannotation__spliceai_max_ds__gte=self.splice_min))
                if self.splice_required and self.splice_allow_null:
                    splicing_q_list.append(Q(variantannotation__spliceai_max_ds__isnull=True))
            else:
                for _, (ds, _) in VariantAnnotation.SPLICEAI_DS_DP.items():
                    splicing_q_list.append(Q(**{f"variantannotation__{ds}__gte": self.splice_min}))
                    if self.splice_required and self.splice_allow_null:
                        splicing_q_list.append(Q(**{f"variantannotation__{ds}__isnull": True}))

        if self.alphagenome_splicing_min is not None and self._filter_column_applied("alphagenome_splicing"):
            splicing_q_list.append(Q(variantannotation__alphagenome_splicing__gte=self.alphagenome_splicing_min))
            if self.splice_required and self.splice_allow_null:
                splicing_q_list.append(Q(variantannotation__alphagenome_splicing__isnull=True))

        if splicing_q_list:
            q_splicing = reduce(operator.or_, splicing_q_list)
            if self.splice_required:
                and_filters.append(q_splicing)
            else:
                or_filters.append(q_splicing)

        if self.columns_version == 1:
            if self.cadd_score_min is not None:
                q_cadd = Q(variantannotation__cadd_phred__gte=self.cadd_score_min)
                if self.cadd_score_required:
                    if self.cadd_score_allow_null:
                        q_cadd |= Q(variantannotation__cadd_phred__isnull=True)
                    and_filters.append(q_cadd)
                else:
                    or_filters.append(q_cadd)

            if self.revel_score_min:
                q_revel = Q(variantannotation__revel_score__gte=self.revel_score_min)
                if self.revel_score_required:
                    if self.revel_score_allow_null:
                        q_revel |= Q(variantannotation__revel_score__isnull=True)
                    and_filters.append(q_revel)
                else:
                    or_filters.append(q_revel)
        elif self.columns_version >= 2:
            # Rank score path predictors
            vav = self.analysis.annotation_version.variant_annotation_version
            for path_prediction in vav.get_rankscore_pathogenic_prediction_funcs():
                path_prediction_min = f"{path_prediction}_min"
                if (score_min := getattr(self, path_prediction_min)) is not None:
                    q_path = Q(**{f"variantannotation__{path_prediction}__gte": score_min})
                    if getattr(self, f"{path_prediction}_required"):
                        if getattr(self, f"{path_prediction}_allow_null"):
                            q_path |= Q(**{f"variantannotation__{path_prediction}__isnull": True})
                        and_filters.append(q_path)
                    else:
                        or_filters.append(q_path)

            # Raw-score filters (v4 onward) — driven by TOOLS so the uncalibrated tools
            # (VARITY_ER, EVE, popEVE, ProtVar ddG, PromoterAI) are offered alongside the
            # calibrated ones.
            if self.columns_version >= 4:
                for tool in TOOLS:
                    raw_field = tool.raw_field
                    if not raw_field or not self._filter_column_applied(raw_field):
                        continue
                    if (threshold := getattr(self, tool.node_threshold_field)) is not None:
                        q_raw = self._raw_score_q(tool, threshold)
                        if getattr(self, f"{raw_field}_required"):
                            if getattr(self, f"{raw_field}_allow_null"):
                                q_raw |= Q(**{f"variantannotation__{raw_field}__isnull": True})
                            and_filters.append(q_raw)
                        else:
                            or_filters.append(q_raw)

                # Categorical pred-field filters
                for pred_field in TOOLS_BY_PRED_FIELD:
                    pred_value = getattr(self, pred_field)
                    if pred_value:
                        q_pred = Q(**{f"variantannotation__{pred_field}": pred_value})
                        if getattr(self, f"{pred_field}_required"):
                            if getattr(self, f"{pred_field}_allow_null"):
                                q_pred |= Q(**{f"variantannotation__{pred_field}__isnull": True})
                            and_filters.append(q_pred)
                        else:
                            or_filters.append(q_pred)

            if self.nmd_escaping_variant:
                q_nmd = Q(variantannotation__nmd_escaping_variant=self.nmd_escaping_variant)
                if self.nmd_escaping_variant_required:
                    and_filters.append(q_nmd)
                else:
                    or_filters.append(q_nmd)

            if self.ptc_nmd and self.columns_version >= 5:
                # Skip on pre-backfill VAVs: nmd_escape_status is null on every row annotated
                # before the PTC calculation existed, so the filter would exclude the whole
                # partition. See VariantAnnotationVersion.backfilled_ptc; get_warnings()
                # surfaces the skip to the user.
                vav = self.analysis.annotation_version.variant_annotation_version
                if vav.backfilled_ptc:
                    nmd_escape_statuses = PTCNMDFilter(self.ptc_nmd).nmd_escape_statuses
                    q_ptc = Q(variantannotation__nmd_escape_status__in=nmd_escape_statuses)
                    if self.ptc_nmd_required:
                        and_filters.append(q_ptc)
                    else:
                        or_filters.append(q_ptc)

            if self.aloft:
                q_aloft = Q(variantannotation__aloft_high_confidence=True)
                ALOFT_OPTIONS = {
                    ALoFTPredictionOptions.DAMAGING: [ALoFTPrediction.RECESSIVE, ALoFTPrediction.DOMINANT],
                }
                q_aloft &= Q(variantannotation__aloft_pred__in=ALOFT_OPTIONS.get(self.aloft, [self.aloft]))
                if self.aloft_required:
                    if self.aloft_allow_null:
                        q_aloft |= Q(variantannotation__aloft_pred__isnull=True)
                    and_filters.append(q_aloft)
                else:
                    or_filters.append(q_aloft)

            if self.splicevault_out_of_frame and self._filter_column_applied("splicevault_out_of_frame_events"):
                top_event_only = self.splicevault_out_of_frame == SpliceVaultOutOfFrame.TOP
                q_splicevault = splicevault_out_of_frame_q(top_event_only)
                if self.splicevault_out_of_frame_required:
                    if self.splicevault_out_of_frame_allow_null:
                        q_splicevault |= Q(variantannotation__splicevault_out_of_frame_events__isnull=True)
                    and_filters.append(q_splicevault)
                else:
                    or_filters.append(q_splicevault)

            if self.utr_high_impact and self._filter_column_applied("utr_annotator_annotation"):
                q_utr = utr_high_impact_q(out_of_frame_only=self.utr_high_impact == UTRHighImpact.OUT_OF_FRAME)
                if self.utr_high_impact_required:
                    if self.utr_high_impact_allow_null:
                        q_utr |= Q(variantannotation__utr_annotator_consequence__isnull=True)
                    and_filters.append(q_utr)
                else:
                    or_filters.append(q_utr)

        if self.damage_predictions_min is not None:
            # Skip the damage_predictions_min filter on pre-backfill VAVs: the
            # predictions_num_pathogenic / predictions_num_benign aggregates default
            # to 0 on legacy rows, so any non-zero damage_predictions_min would
            # silently exclude every row in the partition. See
            # VariantAnnotationVersion.backfilled_damage_counts; get_warnings()
            # surfaces the skip to the user.
            vav = self.analysis.annotation_version.variant_annotation_version
            if vav.backfilled_damage_counts:
                q_damage = Q(variantannotation__predictions_num_pathogenic__gte=self.damage_predictions_min)
                if self.damage_predictions_required:
                    if self.damage_predictions_allow_null:
                        max_benign = self.num_prediction_fields - self.damage_predictions_min
                        q_damage = Q(variantannotation__predictions_num_benign__lte=max_benign)
                    and_filters.append(q_damage)
                else:
                    or_filters.append(q_damage)

        for field, lookup, threshold, required, allow_null in self._gene_constraint_filters():
            q_constraint = gene_constraint_q(field, lookup, threshold, allow_null=required and allow_null)
            if required:
                and_filters.append(q_constraint)
            else:
                or_filters.append(q_constraint)

        if self._applies_outside_repeat_region():
            and_filters.append(Q(variantannotation__repeat_masker__isnull=True))

        for column, threshold, required, allow_null in self._annotsv_threshold_filters():
            q_annotsv = Q(**{f"variantannotation__{column}__gte": threshold})
            if required:
                if allow_null:
                    q_annotsv |= Q(**{f"variantannotation__{column}__isnull": True})
                and_filters.append(q_annotsv)
            else:
                or_filters.append(q_annotsv)

        if self._applies_exclude_sv_artefact_regions():
            and_filters.append(Q(**{f"variantannotation__{column}__isnull": True
                                    for column in ANNOTSV_ARTEFACT_REGION_COLUMNS}))

        if self.cosmic_count_min is not None:
            q_cosmic_count = Q(variantannotation__cosmic_count__gte=self.cosmic_count_min)
            if self.cosmic_count_required:
                and_filters.append(q_cosmic_count)
            else:
                or_filters.append(q_cosmic_count)

        if self.in_cosmic:
            q_in_cosmic = Q(variantannotation__cosmic_id__isnull=False)
            if self.in_cosmic_required:
                and_filters.append(q_in_cosmic)
            else:
                or_filters.append(q_in_cosmic)

        if self.omim_phenotype:
            q_omim = get_omim_q(self.analysis.annotation_version)
            if self.omim_phenotype_required:
                and_filters.append(q_omim)
            else:
                or_filters.append(q_omim)

        if self.protein_domain:
            protein_domains_fields_or = []
            for f in ["interpro_domain", "domains"]:  # TODO: What about UniProt??
                q_domain_field = Q(**{f"variantannotation__{f}__isnull": False})
                protein_domains_fields_or.append(q_domain_field)
            q_protein_domain = reduce(operator.or_, protein_domains_fields_or)
            if self.protein_domain_required:
                and_filters.append(q_protein_domain)
            else:
                or_filters.append(q_protein_domain)

        if self.published:
            published_fields_or = []
            for f in itertools.chain(["pubmed"], VariantAnnotation.MASTERMIND_FIELDS.keys()):
                q_published_field = Q(**{f"variantannotation__{f}__isnull": False})
                published_fields_or.append(q_published_field)

            q_published = reduce(operator.or_, published_fields_or)
            if self.published_required:
                and_filters.append(q_published)
            else:
                or_filters.append(q_published)

        if or_filters:
            q_or = reduce(operator.or_, or_filters)
            and_filters.append(q_or)

        if and_filters:
            q = reduce(operator.and_, and_filters)
        else:
            q = None
        return q

    @property
    def num_prediction_fields(self) -> Optional[int]:
        vav = self.analysis.annotation_version.variant_annotation_version
        if vav is None:
            return None
        if vav.columns_version >= 4:
            return len(vav.get_raw_score_pathogenic_prediction_funcs())
        return len(vav.get_rankscore_pathogenic_prediction_funcs())

    def _get_method_summary(self):
        if self.modifies_parents():
            method_summary = self.get_node_name()
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        """ The headline filters, a line each - everything else is a chip """
        lines = []
        if self.modifies_parents():
            if self.impact_min:
                impact = dict(PathogenicityImpact.CHOICES)[self.impact_min]
                lines.append(f"≥ {impact.capitalize()} impact")
            if self.damage_predictions_min:
                if (num_prediction_fields := self.num_prediction_fields) is None:
                    lines.append(f"{self.damage_predictions_min} damage predictions")
                else:
                    lines.append(f"{self.damage_predictions_min} of {num_prediction_fields} predictions")
            if self.ptc_nmd:
                lines.append(PTCNMDFilter(self.ptc_nmd).label)
        return "\n".join(lines)

    @staticmethod
    def get_help_text() -> str:
        return "Variant type, impact, damage predictions, conservation and splicing filter"

    def get_node_chips(self) -> list[NodeChip]:
        chips = super().get_node_chips()
        if restricts_variant_types(self.variant_types):
            labels = ", ".join(get_variant_type_labels(self.variant_types))
            chips.append(NodeChip(text=labels, icon="fa-solid fa-shapes", title=f"Variant type: {labels}"))
        splice_titles = []
        if self.splice_min is not None:
            splice_titles.append(f"Splicing prediction score >= {self.splice_min}")
        if self.alphagenome_splicing_min is not None:
            splice_titles.append(f"AlphaGenome splicing >= {self.alphagenome_splicing_min}")
        if splice_titles:
            chips.append(NodeChip(text="splice", icon="fa-solid fa-scissors", title="\n".join(splice_titles)))
        if self.splicevault_out_of_frame:
            chips.append(NodeChip(text="SpliceVault", icon="fa-solid fa-scissors",
                                  title=f"SpliceVault: {SpliceVaultOutOfFrame(self.splicevault_out_of_frame).label}"))
        if self.utr_high_impact:
            chips.append(NodeChip(text="5' UTR", title=f"{UTRHighImpact(self.utr_high_impact).label} (UTRAnnotator)"))
        chips += self._score_chips()
        if self.cosmic_count_min is not None:
            chips.append(self._filter_chip("COSMIC", f"COSMIC count ≥ {self.cosmic_count_min}",
                                           self.cosmic_count_required))
        if self.in_cosmic:
            chips.append(self._filter_chip("COSMIC", "Has a COSMIC id", self.in_cosmic_required))
        if self.omim_phenotype:
            chips.append(self._filter_chip("OMIM", "Gene has OMIM phenotypes", self.omim_phenotype_required))
        if self.protein_domain:
            chips.append(self._filter_chip("domain", "In a protein domain", self.protein_domain_required))
        if self.published:
            chips.append(self._filter_chip("published", "Published (PubMed / Mastermind)", self.published_required))
        if self.nmd_escaping_variant:
            chips.append(self._filter_chip("NMD escaping", "NMD escaping variant (LOFTEE)",
                                           self.nmd_escaping_variant_required))
        if self.aloft:
            chips.append(self._filter_chip("ALoFT", f"ALoFT: {ALoFTPredictionOptions(self.aloft).label}",
                                           self.aloft_required, self.aloft_allow_null))
        if self.gene_loeuf_max is not None:
            chips.append(self._filter_chip("LOEUF", f"gnomAD LOEUF ≤ {self.gene_loeuf_max:g}",
                                           self.gene_loeuf_required, self.gene_loeuf_allow_null))
        if self.gene_mis_z_min is not None:
            chips.append(self._filter_chip("missense Z", f"gnomAD missense Z ≥ {self.gene_mis_z_min:g}",
                                           self.gene_mis_z_required, self.gene_mis_z_allow_null))
        if self.outside_repeat_region:
            chips.append(NodeChip(text="no repeat", title="Outside RepeatMasker repeat regions"))
        if self.annotsv_acmg_class_min is not None:
            pathogenicity = Pathogenicity(self.annotsv_acmg_class_min)
            chips.append(self._filter_chip(f"AnnotSV ≥ {pathogenicity.short_label}",
                                           f"AnnotSV ACMG class ≥ {pathogenicity.label}",
                                           self.annotsv_acmg_class_required, self.annotsv_acmg_class_allow_null))
        if self.annotsv_acmg_score_min is not None:
            chips.append(self._filter_chip("AnnotSV score",
                                           f"AnnotSV ranking score ≥ {self.annotsv_acmg_score_min:g}",
                                           self.annotsv_acmg_score_required, self.annotsv_acmg_score_allow_null))
        if self.exclude_sv_artefact_regions:
            chips.append(NodeChip(text="no SV artefact",
                                  title="No SV breakpoint in an ENCODE blacklist, segmental duplication or repeat "
                                        "region (AnnotSV)"))
        return chips

    def _get_node_extra_grid_columns(self) -> list[str]:
        grid_columns = []
        if self._annotsv_threshold_filters():
            grid_columns.append(ANNOTSV_ACMG_GRID_COLUMN)
        if self._applies_exclude_sv_artefact_regions():
            grid_columns.append(ANNOTSV_REGION_GRID_COLUMN)
        return grid_columns

    def _score_chips(self) -> list[NodeChip]:
        """ A chip per tool with a threshold or prediction set - off the TOOLS metadata, so a new tool has one """
        chips = []
        if self.cadd_score_min is not None:  # columns v1 - REVEL's v1 field is the raw tool's
            chips.append(self._filter_chip("CADD", f"CADD phred ≥ {self.cadd_score_min}",
                                           self.cadd_score_required, self.cadd_score_allow_null))
        for tool in TOOLS:
            label = tool.node_label or tool.name
            if tool.raw_field and (value := getattr(self, tool.node_threshold_field)) is not None:
                if tool.raw_direction == RawScoreDirection.LOWER:
                    threshold = f"≤ {value:g}"
                elif tool.raw_direction == RawScoreDirection.MAGNITUDE:
                    threshold = f"|score| ≥ {value:g}"
                else:
                    threshold = f"≥ {value:g}"
                chips.append(self._filter_chip(label, f"{label} {threshold}",
                                               getattr(self, f"{tool.raw_field}_required"),
                                               getattr(self, f"{tool.raw_field}_allow_null")))
            if tool.rankscore_field and (value := getattr(self, f"{tool.rankscore_field}_min")) is not None:
                chips.append(self._filter_chip(f"{tool.name} rank", f"{tool.name} rankscore ≥ {value:g}",
                                               getattr(self, f"{tool.rankscore_field}_required"),
                                               getattr(self, f"{tool.rankscore_field}_allow_null")))
            if tool.pred_field and getattr(self, tool.pred_field):
                prediction = getattr(self, f"get_{tool.pred_field}_display")()
                chips.append(self._filter_chip(tool.name, f"{tool.name}: {prediction}",
                                               getattr(self, f"{tool.pred_field}_required"),
                                               getattr(self, f"{tool.pred_field}_allow_null")))
        return chips

    @staticmethod
    def _filter_chip(text: str, title: str, required: bool, allow_null: Optional[bool] = None) -> NodeChip:
        """ A required filter is marked * - the rest are OR'd together """
        qualifiers = []
        if required:
            text += "*"
            qualifiers.append("required")
        if allow_null:
            qualifiers.append("allows null")
        if qualifiers:
            title += f" ({', '.join(qualifiers)})"
        return NodeChip(text=text, title=title)

    @staticmethod
    def get_node_class_label():
        return "Effect"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-effect")


auditlog.register(DamageNode)
