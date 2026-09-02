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
from snpdb.models import CustomColumn, CustomColumnsCollection, Variant
from snpdb.models.models_enums import ColumnAnnotationLevel
from snpdb.views.datatable_view import CellData, FilterField, NullOrder, RichColumn

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


def _catalogue_column_kwargs(column) -> dict:
    """ RichColumn kwargs a VariantGridColumn carries whether it's shown or riding along hidden - a
        hidden member still gets a CSV header, and the catalogue label beats the model verbose_name
        (locus__ref__seq and alt__seq are both 'seq') """
    kwargs = {
        "model_field": column.model_field,
        "queryset_field": column.queryset_field or column.variant_column in VARIANT_GRID_EXTRA_ANNOTATION_ALIASES,
        "label": column.label,
        "header_title": column.description.replace("'", "&#146;"),
    }
    if column.width is not None:
        kwargs["width"] = column.width
    return kwargs


def _composite_column_kwargs(column, members: list, column_overrides: dict[str, dict]) -> dict:
    """ The extras a composite cell needs on top of its catalogue kwargs - the members it draws, the
        sort keys its header offers, and the generic renderer unless an override names a bespoke one.
        `members` are the ones annotated for this version, so nothing named here is missing from the row """
    def member_entry(member) -> dict:
        entry = {"path": member.column.variant_column, "label": member.column.label}
        if renderer := column_overrides.get(member.column.variant_column, {}).get("client_renderer"):
            entry["renderer"] = renderer
        return entry

    kwargs = {
        "client_renderer": "VariantGridFormat.composite",
        "client_renderer_kwargs": {"members": [member_entry(m) for m in members]},
        "sort_menu": [{"label": m.column.label, "column": m.column.variant_column}
                      for m in members if m.in_sort_menu],
        "orderable": True,
    }
    if not column.queryset_field:
        # Display only - nothing of its own to sort on, so it sorts on the value the cell reads as
        kwargs["sort_keys"] = [members[0].column.variant_column]
    return kwargs


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
    columns_queryset = columns_queryset.prefetch_related("column__composite_members__column")

    def annotated_for_this_version(column) -> bool:
        return column.pk not in ever_referenced or column.pk in in_version_vgcs

    shown_columns = []  # (column, the members it draws that this version annotates)
    for c in columns_queryset:
        column = c.column
        if not column.variant_column:
            continue
        # Tags are only shown in the analysis they are in (otherwise will just show tags_global)
        if column.variant_column == "tags" and not analysis_tags:
            continue
        members = [m for m in column.composite_members.all() if annotated_for_this_version(m.column)]
        if column.is_composite and not members:
            continue  # nothing this version annotates for the cell to draw
        shown_columns.append((column, members))
    shown_paths = {column.variant_column for column, _ in shown_columns}

    fields_kwargs: dict[str, dict] = {}  # field path -> RichColumn kwargs, in column order
    sample_columns_position = None

    for field_pos, (column, members) in enumerate(shown_columns):
        if column.model_field is False:
            if column.annotation_level == ColumnAnnotationLevel.SAMPLE_LEVEL:
                sample_columns_position = field_pos

        kwargs = _catalogue_column_kwargs(column)
        if members:
            kwargs.update(_composite_column_kwargs(column, members, column_overrides))
        fields_kwargs[column.variant_column] = kwargs

        # The members the cell draws ride along hidden, right after it - unless the collection shows
        # one standalone, where it is simply visible, once, in its own place
        for member in members:
            if member.column.variant_column in shown_paths:
                continue
            fields_kwargs[member.column.variant_column] = {**_catalogue_column_kwargs(member.column),
                                                           "visible": False}

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
