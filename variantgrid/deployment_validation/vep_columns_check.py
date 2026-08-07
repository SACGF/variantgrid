def check_vep_columns_registry() -> dict:
    """ Validates annotation/vep_columns.py against VariantGridColumn + internal uniqueness. """
    from annotation.vep_columns import VEP_COLUMNS, all_variant_grid_column_ids
    from snpdb.models import VariantGridColumn

    known = set(VariantGridColumn.objects.values_list("pk", flat=True))
    referenced = all_variant_grid_column_ids()
    missing = sorted(referenced - known)

    seen = set()
    duplicates = []
    for c in VEP_COLUMNS:
        k = (c.source_field, c.vep_plugin, c.vep_custom,
             c.min_columns_version, c.max_columns_version,
             c.min_vep_version, c.max_vep_version,
             c.genome_builds, c.pipeline_types,
             c.variant_grid_columns)
        if k in seen:
            duplicates.append(k)
        seen.add(k)

    data = {
        "unknown_variant_grid_columns": {
            "valid": not missing,
            "fix": f"vep_columns.py references unknown VariantGridColumn ids: {missing}" if missing else None,
        },
        "duplicate_vep_column_defs": {
            "valid": not duplicates,
            "fix": f"Duplicate VEPColumnDef keys: {duplicates}" if duplicates else None,
        },
    }
    return data
