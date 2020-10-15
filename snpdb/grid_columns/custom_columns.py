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
