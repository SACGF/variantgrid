import logging

from django.core.exceptions import FieldError

from snpdb.models import Variant, VariantGridColumn


def check_variantgrid_columns() -> dict:
    # These are added via annotate so won't be in a standard variant query
    # TODO: could build proper querby via get_variantgrid_extra_annotate
    SPECIAL_COLUMNS = {
        "",  # Blank
        "internally_classified", "internally_classified_labs", "max_internal_classification",
        "tags", "tags_global",
        "global_variant_zygosity__hom_count", "global_variant_zygosity__ref_count",
        "global_variant_zygosity__unk_count", "global_variant_zygosity__het_count",
    }

    good_columns = []
    bad_columns = []
    qs = Variant.objects.filter(pk=1)  # So it's fast
    vgc_qs = VariantGridColumn.objects.all()
    vgc_qs = vgc_qs.exclude(variant_column__in=SPECIAL_COLUMNS)
    for vc in vgc_qs.values_list("variant_column", flat=True):
        try:
            qs.filter(**{f"{vc}__isnull": False})
            good_columns.append(vc)
        except FieldError as e:
            logging.error(e)
            bad_columns.append(vc)

    if bad_columns:
        cols = ", ".join(bad_columns)
        fix = f"The following VariantGridColumns.column_name are not valid queryset paths from Variant: '{cols}'"
    else:
        fix = None

    data = {
        "variant_grid_columns": {
            "valid": not bad_columns,
            "fix": fix,
        }
    }
    return data
