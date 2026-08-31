from django.contrib.auth.models import User
from django.db.models import F, Max, OuterRef, Q, StringAgg, Subquery, TextField, Value
from django.db.models.fields.json import KeyTextTransform, KeyTransform
from django.db.models.functions import Cast, Coalesce, Concat, TruncDate

from analysis.models import VariantTag
from annotation import vep_columns
from annotation.models import AnnotationVersion
from classification.models import ClassificationModification
from library.jqgrid.jqgrid_sql import get_overrides
from snpdb.models import CustomColumn, CustomColumnsCollection, VariantGridColumn
from snpdb.models.models_enums import ColumnAnnotationLevel

# Row fields the client renderers read whichever columns the collection shows - they ride along hidden.
# @see VariantGridFormat.representativeVariant / VariantGridFormat.classifications in variantgrid_formats.js
VARIANT_COLUMN_ROW_FIELDS = [
    "locus__contig__name",
    "locus__position",
    "locus__ref__seq",
    "alt__seq",
    "svlen",
    "variantannotation__symbol",
    "variantannotation__hgvs_c",
    "variantannotation__hgvs_p",
    "variantannotation__hgvs_g",
]
CLASSIFICATIONS_COLUMN_ROW_FIELDS = [
    "clinvar__highest_pathogenicity",
    "clinvar__clinical_significance",
    "clinvar__review_status",
    "clinvar__highest_oncogenicity",
    "clinvar__oncogenic_classification",
    "clinvar__oncogenic_review_status",
    "clinvar__somatic_tier",
    "clinvar__somatic_review_status",
]
# Composite cells: the visible column the renderer sits on -> the partner values it draws alongside.
# The partners ride along hidden, and only while the collection actually shows the column that reads
# them. @see _get_standard_overrides in snpdb/grids.py for the renderer each one is paired with
COMPOSITE_COLUMN_ROW_FIELDS = {
    "variantannotation__consequence": ["variantannotation__impact"],
    "variantannotation__gnomad_af": ["variantannotation__gnomad_filtered"],
    "variantannotation__gnomad_popmax_af": ["variantannotation__gnomad_popmax"],
    "variantannotation__spliceai_max_ds": [
        "variantannotation__spliceai_pred_ds_ag",
        "variantannotation__spliceai_pred_ds_al",
        "variantannotation__spliceai_pred_ds_dg",
        "variantannotation__spliceai_pred_ds_dl",
        "variantannotation__spliceai_pred_dp_ag",
        "variantannotation__spliceai_pred_dp_al",
        "variantannotation__spliceai_pred_dp_dg",
        "variantannotation__spliceai_pred_dp_dl",
    ],
    "variantannotation__maxentscan_percent_diff_ref": [
        "variantannotation__maxentscan_ref",
        "variantannotation__maxentscan_alt",
        "variantannotation__maxentscan_diff",
    ],
    "variantannotation__mastermind_count_1_cdna": [
        "variantannotation__mastermind_count_2_cdna_prot",
        "variantannotation__mastermind_count_3_aa_change",
        "variantannotation__mastermind_mmid3",
    ],
    "variantannotation__aloft_pred": [
        "variantannotation__aloft_high_confidence",
        "variantannotation__aloft_prob_tolerant",
        "variantannotation__aloft_prob_recessive",
        "variantannotation__aloft_prob_dominant",
        "variantannotation__aloft_ensembl_transcript",
    ],
    "variantannotation__predictions_num_pathogenic": ["variantannotation__predictions_num_benign"],
    "global_variant_zygosity__het_count": [
        "global_variant_zygosity__hom_count",
        "global_variant_zygosity__ref_count",
        "global_variant_zygosity__unk_count",
    ],
}
# These come from get_variantgrid_extra_annotate rather than the model, so they are colmodel-only
# (not in the values() queryset)
CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS = [
    "internally_classified",
    "max_internal_classification",
    "internally_classified_somatic",
    "max_internal_somatic_classification",
]


def get_custom_column_fields_override_and_sample_position(custom_columns_collection: CustomColumnsCollection,
                                                          annotation_version: AnnotationVersion,
                                                          analysis_tags=False):
    variant_annotation_version = annotation_version.variant_annotation_version
    # A column only annotated for other builds (eg the gnomAD 4 FAFs, GRCh38 only) would always be
    # empty here, so it drops out of the grid (and its exports) for this build
    in_version_vgcs = {
        vgc for c in vep_columns.filter_for(columns_version=variant_annotation_version.columns_version,
                                            genome_build_name=variant_annotation_version.genome_build.name)
        for vgc in c.variant_grid_columns
    }
    ever_referenced = vep_columns.all_variant_grid_column_ids()
    q_columns_this_version = ~Q(column__in=ever_referenced) | Q(column__in=in_version_vgcs)
    columns_queryset = CustomColumn.objects.filter(q_columns_this_version,
                                                   custom_columns_collection=custom_columns_collection)
    columns_queryset = columns_queryset.select_related("column").order_by("sort_order").distinct()
    fields = []
    sample_columns_position = None
    override = {}

    for field_pos, c in enumerate(columns_queryset):
        if f := c.column.variant_column:
            # Tags are only shown in the analysis they are in (otherwise will just show tags_global)
            if f == "tags" and not analysis_tags:
                continue
            fields.append(f)

        if c.column.model_field is False:
            if c.column.annotation_level == ColumnAnnotationLevel.SAMPLE_LEVEL:
                sample_columns_position = field_pos

        col_overrides = get_overrides([f], [{}],
                                      model_field=c.column.model_field, queryset_field=c.column.queryset_field)
        col_override = col_overrides[f]
        description = c.column.description.replace("'", "&#146;")
        col_override["headerTitle"] = description
        col_override["label"] = c.column.label

        if c.column.width is not None:
            col_override["width"] = c.column.width

        override[f] = col_override

    composite_partners = [partner for visible, partners in COMPOSITE_COLUMN_ROW_FIELDS.items()
                          for partner in partners if visible in fields]
    hidden_fields = VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS + composite_partners
    # A hidden field still gets a CSV header - use the catalogue label rather than the model
    # verbose_name (locus__ref__seq and alt__seq are both 'seq')
    hidden_columns = {vgc.variant_column: vgc
                      for vgc in VariantGridColumn.objects.filter(variant_column__in=hidden_fields)}
    for field in hidden_fields:
        if field not in fields:
            fields.append(field)
            vgc = hidden_columns.get(field)
            if vgc is not None and vgc.model_field is False:
                # An annotation rather than a model field (the global zygosity counts) - the colmodel
                # has to say so, or the engine goes looking for a model field of that name
                ov = get_overrides([field], [{}], model_field=False,
                                   queryset_field=vgc.queryset_field)[field]
            else:
                ov = {}
            if vgc is not None:
                ov["label"] = vgc.label
            ov["hidden"] = True
            override[field] = ov

    # Grid-only values from get_variantgrid_extra_annotate - not queryset fields, so they need a
    # colmodel that says so
    for field in CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS:
        if field not in fields:
            fields.append(field)
            ov = get_overrides([field], [{}], model_field=False, queryset_field=False)[field]
            ov["hidden"] = True
            override[field] = ov

    return fields, override, sample_columns_position


def get_variantgrid_extra_annotate(user: User, exclude_analysis=None) -> dict:
    # Ordering must be cleared before the .values() grouping - ordered columns get added to GROUP BY,
    # which would split the per-allele aggregates into one group per record. The StringAgg order_by
    # keeps the parallel '|' separated fields lining up
    classification_qs = ClassificationModification.latest_for_user(user).filter(classification__allele__variantallele__variant_id=OuterRef("id")).order_by()
    by_allele = classification_qs.values("classification__allele")
    internally_classified = by_allele.annotate(cs_summary=StringAgg(Coalesce("classification__clinical_significance", Value('U')), delimiter=Value('|'), order_by="pk")).values_list("cs_summary")
    internally_classified_labs = by_allele.annotate(c_lab=StringAgg(Coalesce("classification__lab__name", Value('')), delimiter=Value('|'), order_by="pk")).values_list("c_lab")
    max_internal_classification = by_allele.annotate(cs_max=Max("classification__clinical_significance")).values_list("cs_max")

    somatic_cs = KeyTextTransform("clinical_significance", KeyTransform("somatic", "classification__summary"))
    internally_classified_somatic = by_allele.annotate(scs_summary=StringAgg(Coalesce(somatic_cs, Value('U'), output_field=TextField()), delimiter=Value('|'), order_by="pk")).values_list("scs_summary")
    # Somatic tiers don't sort lexically (tier_1_or_2 between tier_1/tier_2) so take the value of the
    # highest summary sort score rather than Max
    max_internal_somatic_classification = classification_qs.annotate(scs=somatic_cs).exclude(scs__isnull=True).order_by(F("classification__summary__somatic__sort").desc(nulls_last=True)).values_list("scs")

    tags_qs = VariantTag.filter_for_user(user).filter(allele__variantallele__variant_id=OuterRef("id"))
    if exclude_analysis:
        tags_qs = tags_qs.filter(Q(analysis__isnull=True) | Q(analysis__id__ne=exclude_analysis.pk))
    # "tag_id:date" entries, e.g. "Artefact:2024-03-01|SomaticReportable:2023-05-06" - the date lets
    # formatters show fresh vs total counts (see tagsGlobalFormatter / format_items_iterator)
    tag_and_date = Concat("tag_id", Value(":"), Cast(TruncDate("created"), TextField()), output_field=TextField())
    tags_global = tags_qs.values("allele").annotate(tags=StringAgg(tag_and_date, delimiter=Value('|'))).values_list("tags")

    return {
        "internally_classified": Subquery(internally_classified[:1]),
        "internally_classified_labs": Subquery(internally_classified_labs[:1]),
        "max_internal_classification": Subquery(max_internal_classification[:1]),
        "internally_classified_somatic": Subquery(internally_classified_somatic[:1]),
        "max_internal_somatic_classification": Subquery(max_internal_somatic_classification[:1]),
        "tags_global": Subquery(tags_global[:1]),
    }
