from django.contrib.auth.models import User

from library.jqgrid_sql import get_overrides
from snpdb.models import CustomColumn, CustomColumnsCollection
from snpdb.models.models_enums import ColumnAnnotationLevel


def get_custom_column_fields_override_and_sample_position(custom_columns_collection: CustomColumnsCollection):
    columns_queryset = CustomColumn.objects.filter(custom_columns_collection=custom_columns_collection).order_by("sort_order")
    fields = []
    sample_columns_position = None
    override = {}

    for field_pos, c in enumerate(columns_queryset):
        if f := c.column.variant_column:
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

    # Used by detailsLink() in JavaScript
    ID_FORMATTER_REQUIRED_FIELDS = ["locus__contig__name",
                                    "locus__position",
                                    "clinvar__highest_pathogenicity",
                                    "clinvar__clinical_significance"]
    for field in ID_FORMATTER_REQUIRED_FIELDS:
        if not field in fields:
            fields.append(field)
            ov = override.get(field, {})
            ov["hidden"] = True
            override[field] = ov

    return fields, override, sample_columns_position


SELECT_INTERNALLY_CLASSIFIED_SQL = """
select string_agg(coalesce(classification_classification.clinical_significance, 'U'), '|')
from classification_classification
where classification_classification.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where
    allele_id in (
        select allele_id from snpdb_variantallele where variant_id = snpdb_variant.id
    )
)
"""

SELECT_MAX_INTERNAL_CLASSIFICATION = """
select max(coalesce(classification_classification.clinical_significance, '0'))
from classification_classification
where classification_classification.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where
    allele_id in (
        select allele_id from snpdb_variantallele where variant_id = snpdb_variant.id
    )
)
"""


SELECT_TAGGED_SQL = """
select string_agg(coalesce(analysis_varianttag.tag_id, 'U'), '|')
from analysis_varianttag
where analysis_varianttag.variant_id in (
    select snpdb_variantallele.variant_id
    from snpdb_variantallele
    where
    allele_id in (
        select allele_id from snpdb_variantallele where variant_id = snpdb_variant.id
    )
)
"""


def get_variantgrid_extra_alias_and_select_columns(user: User, exclude_analysis=None):
    # TODO: Need to add user level security to classifications and tags
    tags_global_sql = SELECT_TAGGED_SQL
    if exclude_analysis:
        tags_global_sql += " AND analysis_varianttag.analysis_id <> %d" % exclude_analysis.pk

    alias_and_select = {
        "internally_classified": SELECT_INTERNALLY_CLASSIFIED_SQL,
        "max_internal_classification": SELECT_MAX_INTERNAL_CLASSIFICATION,
        "tags_global": tags_global_sql,
    }
    return alias_and_select.items()
