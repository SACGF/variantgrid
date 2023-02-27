from typing import Dict

from django.contrib.auth.models import User
from django.contrib.postgres.aggregates import StringAgg
from django.db.models import Q, OuterRef, Subquery, Max, Value
from django.db.models.functions import Coalesce

from analysis.models import VariantTag
from annotation.models import AnnotationVersion, ColumnVEPField
from classification.models import Classification, ClassificationModification
from library.jqgrid.jqgrid_sql import get_overrides
from snpdb.models import CustomColumn, CustomColumnsCollection
from snpdb.models.models_enums import ColumnAnnotationLevel


def get_custom_column_fields_override_and_sample_position(custom_columns_collection: CustomColumnsCollection,
                                                          annotation_version: AnnotationVersion):
    q_cvf = ColumnVEPField.get_columns_version_q(annotation_version.variant_annotation_version.columns_version)
    cvf_qs = ColumnVEPField.objects.filter(q_cvf)
    q_columns_this_version = Q(column__columnvepfield__isnull=True) | Q(column__columnvepfield__in=cvf_qs)
    columns_queryset = CustomColumn.objects.filter(q_columns_this_version,
                                                   custom_columns_collection=custom_columns_collection)
    columns_queryset = columns_queryset.select_related("column").order_by("sort_order").distinct()
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


def get_variantgrid_extra_annotate(user: User, exclude_analysis=None) -> Dict:
    classification_qs = ClassificationModification.latest_for_user(user).filter(classification__allele__variantallele__variant_id=OuterRef("id"))
    internally_classified = classification_qs.annotate(cs=Coalesce("classification__clinical_significance", Value('U'))).values("classification__allele").annotate(cs_summary=StringAgg("cs", delimiter='|')).values_list("cs_summary")
    max_internal_classification = classification_qs.annotate(cs=Coalesce("classification__clinical_significance", Value('0'))).values("classification__allele").annotate(cs_max=Max("classification__clinical_significance")).values_list("cs_max")

    tags_qs = VariantTag.filter_for_user(user).filter(allele__variantallele__variant_id=OuterRef("id"))
    if exclude_analysis:
        tags_qs = tags_qs.filter(Q(analysis__isnull=True) | Q(analysis__id__ne=exclude_analysis.pk))
    tags_global = tags_qs.values("allele").annotate(tags=StringAgg("tag_id", delimiter='|')).values_list("tags")

    return {
        "internally_classified": Subquery(internally_classified[:1]),
        "max_internal_classification": Subquery(max_internal_classification[:1]),
        "tags_global": Subquery(tags_global[:1]),
    }
