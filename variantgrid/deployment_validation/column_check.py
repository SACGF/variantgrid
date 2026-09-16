import logging

from django.core.exceptions import FieldError

from snpdb.models import Variant, VariantGridColumn

# Display only columns that are neither a composite nor drawn inside one - tags come from the
# analysis, Sample marks where the per-sample columns go
_DISPLAY_ONLY_WITHOUT_MEMBERS = {"tags", "tags_global", "Sample"}


def check_variantgrid_columns() -> dict:
    # The zygosity counts come from VariantZygosityCountCollection.annotate, so they aren't a path
    # on a standard variant query either - every other annotate alias is queryset_field=False
    # TODO: could build proper querby via get_variantgrid_extra_annotate
    SPECIAL_COLUMNS = {
        "global_variant_zygosity__hom_count", "global_variant_zygosity__ref_count",
        "global_variant_zygosity__unk_count", "global_variant_zygosity__het_count",
    }

    good_columns = []
    bad_columns = []
    qs = Variant.objects.filter(pk=1)  # So it's fast
    # A display only column (a composite, tags) draws its value from a renderer - it has no path
    vgc_qs = VariantGridColumn.objects.exclude(queryset_field=False)
    vgc_qs = vgc_qs.exclude(variant_column__in=SPECIAL_COLUMNS)
    for vc in vgc_qs.values_list("variant_column", flat=True):
        try:
            qs.filter(**{f"{vc}__isnull": False})
            good_columns.append(vc)
        except FieldError as e:
            logging.error(e)
            bad_columns.append(vc)

    # A display only column has to be a composite - it has nothing else to draw. The exceptions are
    # the annotate aliases a composite draws (they're members) and the ones with their own renderer
    memberless = list(VariantGridColumn.objects.filter(queryset_field=False, composite_members__isnull=True,
                                                       composite_membership__isnull=True)
                      .exclude(pk__in=_DISPLAY_ONLY_WITHOUT_MEMBERS).values_list("pk", flat=True))

    fixes = []
    if bad_columns:
        cols = ", ".join(bad_columns)
        fixes.append(f"The following VariantGridColumns.column_name are not valid queryset paths from "
                     f"Variant: '{cols}'")
    if memberless:
        cols = ", ".join(memberless)
        fixes.append(f"The following display only VariantGridColumns have no CompositeColumnMember "
                     f"rows, so they would draw nothing: '{cols}'")

    data = {
        "variant_grid_columns": {
            "valid": not fixes,
            "fix": "\n".join(fixes) or None,
        }
    }
    return data
