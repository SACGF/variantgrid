"""
Maps AnnotSV TSV output fields -> VariantGridColumn destinations.

The TSV inserter uses this to know what to write, and the annotation descriptions
page uses it to show which columns came from AnnotSV rather than VEP.
"""
from dataclasses import dataclass
from typing import Optional

from annotation.models.models_enums import ColumnAnnotationCategory


@dataclass(frozen=True)
class AnnotSVColumnDef:
    """ `tsv_column` False = assembled from several TSV columns by the inserter rather than read directly """
    source_field: str
    variant_grid_column: str
    category: ColumnAnnotationCategory
    source_field_processing_description: Optional[str] = None
    tsv_column: bool = True


ANNOTSV_COLUMNS: tuple[AnnotSVColumnDef, ...] = (
    AnnotSVColumnDef("ACMG_class", "annotsv_acmg_class", ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS),
    AnnotSVColumnDef("AnnotSV_ranking_score", "annotsv_acmg_score",
                     ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS),
    AnnotSVColumnDef("AnnotSV_ranking_criteria", "annotsv_acmg_criteria",
                     ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS),
    AnnotSVColumnDef("RE_gene", "annotsv_re_gene", ColumnAnnotationCategory.GENE_ANNOTATIONS),
    AnnotSVColumnDef("Repeat_type_left", "annotsv_repeat_type_left", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("Repeat_type_right", "annotsv_repeat_type_right", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("SegDup_left", "annotsv_segdup_left", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("SegDup_right", "annotsv_segdup_right", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("ENCODE_blacklist_left", "annotsv_encode_blacklist_left", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("ENCODE_blacklist_right", "annotsv_encode_blacklist_right", ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("ENCODE_blacklist_characteristics_left", "annotsv_encode_blacklist_characteristics_left",
                     ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("ENCODE_blacklist_characteristics_right", "annotsv_encode_blacklist_characteristics_right",
                     ColumnAnnotationCategory.SEQUENCE),
    AnnotSVColumnDef("B_gain_AFmax", "annotsv_b_gain_af_max", ColumnAnnotationCategory.FREQUENCY_DATA),
    AnnotSVColumnDef("B_loss_AFmax", "annotsv_b_loss_af_max", ColumnAnnotationCategory.FREQUENCY_DATA),
    AnnotSVColumnDef("B_ins_AFmax", "annotsv_b_ins_af_max", ColumnAnnotationCategory.FREQUENCY_DATA),
    AnnotSVColumnDef("B_inv_AFmax", "annotsv_b_inv_af_max", ColumnAnnotationCategory.FREQUENCY_DATA),
    AnnotSVColumnDef("Frameshift", "annotsv_frameshift", ColumnAnnotationCategory.FUNCTIONAL_EFFECT),
    AnnotSVColumnDef("Exons_spanned", "annotsv_exons_spanned", ColumnAnnotationCategory.FUNCTIONAL_EFFECT),
    AnnotSVColumnDef("Dist_nearest_SS", "annotsv_dist_nearest_ss", ColumnAnnotationCategory.SPLICING_PREDICTIONS),
    AnnotSVColumnDef("Nearest_SS_type", "annotsv_nearest_ss_type", ColumnAnnotationCategory.SPLICING_PREDICTIONS),
    AnnotSVColumnDef("OMIM_inheritance", "annotsv_omim_inheritance", ColumnAnnotationCategory.PHENOTYPE),
    AnnotSVColumnDef("OMIM_morbid", "annotsv_omim_morbid", ColumnAnnotationCategory.PHENOTYPE),
    AnnotSVColumnDef("OMIM_phenotype", "annotsv_omim_phenotype", ColumnAnnotationCategory.PHENOTYPE),
    AnnotSVColumnDef("OMIM_ID", "annotsv_omim_id", ColumnAnnotationCategory.PHENOTYPE),
    AnnotSVColumnDef("P_{gain,loss,ins,inv}_{source,phen,hpo,coord}", "annotsv_pathogenic_overlaps",
                     ColumnAnnotationCategory.PATHOGENICITY_PREDICTIONS,
                     source_field_processing_description="collapsed into JSON, keyed by event type",
                     tsv_column=False),
)

# AnnotSV TSV column name -> VariantAnnotation field name, for the fields read straight off a row.
# Values verified against the bundle shipped with AnnotSV 3.5.8.
FULL_COLUMN_MAP: dict[str, str] = {
    c.source_field: c.variant_grid_column for c in ANNOTSV_COLUMNS if c.tsv_column
}


def for_variant_grid_column(vgc_id: str) -> Optional[AnnotSVColumnDef]:
    for c in ANNOTSV_COLUMNS:
        if c.variant_grid_column == vgc_id:
            return c
    return None


def all_variant_grid_column_ids() -> frozenset[str]:
    return frozenset(c.variant_grid_column for c in ANNOTSV_COLUMNS)
