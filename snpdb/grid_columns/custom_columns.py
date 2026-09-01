from typing import Optional

from django.contrib.auth.models import User
from django.db.models import F, Max, OuterRef, Q, StringAgg, Subquery, TextField, Value, fields
from django.db.models.fields.json import KeyTextTransform, KeyTransform
from django.db.models.functions import Cast, Coalesce, Concat, TruncDate
from django.utils.timezone import localtime

from analysis.models import VariantTag
from annotation import vep_columns
from annotation.models import AnnotationVersion
from classification.models import ClassificationModification
from library.django_utils import resolve_field_path
from library.utils import JsonDataType, add_exception_note
from snpdb.models import CustomColumn, CustomColumnsCollection, Variant, VariantGridColumn
from snpdb.models.models_enums import ColumnAnnotationLevel
from snpdb.views.datatable_view import CellData, FilterField, NullOrder, RichColumn

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
# These come from get_variantgrid_extra_annotate rather than the model
CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS = [
    "internally_classified",
    "max_internal_classification",
    "internally_classified_somatic",
    "max_internal_somatic_classification",
]
# Every alias get_variantgrid_extra_annotate adds - they're selected out like a model field would be,
# they just have no Django field to introspect for a label or a filter
VARIANT_GRID_EXTRA_ANNOTATION_ALIASES = set(CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS) | {
    "internally_classified_labs",
    "tags_global",
}


# Filter type offered for a model field's class (first match wins)
_FIELD_FILTER_TYPES = [
    (fields.AutoField, 'int'),
    (fields.IntegerField, 'int'),
    (fields.FloatField, 'float'),
    (fields.DateTimeField, 'date'),
]


def _make_choices_renderer(choices: dict):
    """ Need this to perform a closure over a loop variable """
    def choices_renderer(cell: CellData) -> JsonDataType:
        value = cell.value
        return choices.get(value, value)
    return choices_renderer


def _render_local_datetime(cell: CellData) -> JsonDataType:
    if value := cell.value:
        return localtime(value).strftime("%Y-%m-%d %H:%M")
    return value


def _model_field_column_kwargs(field) -> dict:
    """ RichColumn kwargs derived from the Django field a column path resolves to - its label, the
        filter the client should offer, and the server side formatting the CSV shares with the grid """
    kwargs = {"label": field.verbose_name}
    if isinstance(field, fields.related.ForeignKey):
        # A rule on the FK itself would be a type error - use the full path to the final column
        return kwargs

    if isinstance(field, fields.BooleanField):
        kwargs["column_filter"] = FilterField('select', {'False': 'False', 'True': 'True'})
    elif field.choices:
        kwargs["column_filter"] = FilterField('select', dict(field.choices))
    else:
        filter_type = 'text'
        for field_class, name in _FIELD_FILTER_TYPES:
            if isinstance(field, field_class):
                filter_type = name
                break
        kwargs["column_filter"] = FilterField(filter_type)

    if isinstance(field, fields.CharField) and field.choices:
        kwargs["renderer"] = _make_choices_renderer(dict(field.choices))
        kwargs["csv_rendered"] = True
    elif isinstance(field, fields.DateTimeField):
        kwargs["renderer"] = _render_local_datetime
        kwargs["csv_rendered"] = True
    return kwargs


def variant_column_rich_column(path: str, model_field: bool = True, queryset_field: bool = True,
                               **overrides) -> RichColumn:
    """ A variant grid column for a `VariantGridColumn.variant_column` path.

        model_field: resolve the path through the model for its label, filter and formatting.
        queryset_field: the path (or annotation alias) is selected out - display only columns take
        their value from a renderer instead. """
    kwargs = {
        "key": path if queryset_field else None,
        "name": path,
        "orderable": True,
        "null_order": NullOrder.FIRST_ON_ASC,
        "search": False,  # the variant grids have no search box - filtering is the filter builder's
        "include_in_csv": True,
    }
    if model_field:
        kwargs.update(_model_field_column_kwargs(resolve_field_path(Variant._meta, path)))
    kwargs.update(overrides)
    if not kwargs.get("key") and not kwargs.get("sort_keys"):
        kwargs["orderable"] = False  # nothing to order on
    try:
        return RichColumn(**kwargs)
    except ValueError as ve:
        add_exception_note(ve, f"Building variant grid column '{path}'")
        raise


def get_variant_grid_columns(custom_columns_collection: CustomColumnsCollection,
                             annotation_version: AnnotationVersion,
                             column_overrides: dict[str, dict],
                             analysis_tags=False) -> tuple[list[RichColumn], Optional[int]]:
    """ (the collection's columns plus the hidden ones the renderers read, where sample columns go) """
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

    fields_kwargs: dict[str, dict] = {}  # field path -> RichColumn kwargs, in column order
    sample_columns_position = None

    for field_pos, c in enumerate(columns_queryset):
        f = c.column.variant_column
        if not f:
            continue
        # Tags are only shown in the analysis they are in (otherwise will just show tags_global)
        if f == "tags" and not analysis_tags:
            continue

        if c.column.model_field is False:
            if c.column.annotation_level == ColumnAnnotationLevel.SAMPLE_LEVEL:
                sample_columns_position = field_pos

        fields_kwargs[f] = {
            "model_field": c.column.model_field,
            "queryset_field": c.column.queryset_field or f in VARIANT_GRID_EXTRA_ANNOTATION_ALIASES,
            "label": c.column.label,
            "header_title": c.column.description.replace("'", "&#146;"),
        }
        if c.column.width is not None:
            fields_kwargs[f]["width"] = c.column.width

    composite_partners = [partner for visible, partners in COMPOSITE_COLUMN_ROW_FIELDS.items()
                          for partner in partners if visible in fields_kwargs]
    hidden_fields = VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS + composite_partners
    # A hidden field still gets a CSV header - use the catalogue label rather than the model
    # verbose_name (locus__ref__seq and alt__seq are both 'seq')
    hidden_columns = {vgc.variant_column: vgc
                      for vgc in VariantGridColumn.objects.filter(variant_column__in=hidden_fields)}
    for path in hidden_fields:
        if path in fields_kwargs:
            continue
        kwargs = {"visible": False}
        if vgc := hidden_columns.get(path):
            kwargs["label"] = vgc.label
            kwargs["model_field"] = vgc.model_field
            kwargs["queryset_field"] = vgc.queryset_field
        fields_kwargs[path] = kwargs

    # Grid only values from get_variantgrid_extra_annotate - annotation aliases, not model fields
    for path in CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS:
        if path not in fields_kwargs:
            fields_kwargs[path] = {"visible": False, "model_field": False}

    rich_columns = []
    for path, kwargs in fields_kwargs.items():
        kwargs = {**kwargs, **column_overrides.get(path, {})}
        rich_columns.append(variant_column_rich_column(path, **kwargs))
    return rich_columns, sample_columns_position


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
