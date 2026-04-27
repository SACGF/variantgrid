from django.core.checks import Error, Tags, register


@register(Tags.models)
def check_vep_columns_registry(app_configs, **kwargs):
    """ Validates annotation/vep_columns.py against VariantGridColumn + internal uniqueness. """
    from annotation.vep_columns import VEP_COLUMNS, all_variant_grid_column_ids
    from snpdb.models import VariantGridColumn

    errors = []

    try:
        known = set(VariantGridColumn.objects.values_list("pk", flat=True))
    except Exception:
        # DB not ready (e.g. during initial migrate); skip — system check will run again later.
        return []

    referenced = all_variant_grid_column_ids()
    if missing := referenced - known:
        errors.append(Error(
            f"vep_columns.py references unknown VariantGridColumn ids: {sorted(missing)}",
            id="annotation.E001",
        ))

    seen = set()
    for c in VEP_COLUMNS:
        k = (c.source_field, c.vep_plugin, c.vep_custom,
             c.min_columns_version, c.max_columns_version,
             c.min_vep_version, c.max_vep_version,
             c.genome_builds, c.pipeline_types,
             c.variant_grid_columns)
        if k in seen:
            errors.append(Error(
                f"Duplicate VEPColumnDef key: {k}",
                id="annotation.E002",
            ))
        seen.add(k)

    return errors
