
def get_overrides(columns, column_data, model_field=True, queryset_field=True):
    base_colmodel_override = {
        'model_field': True,
        'queryset_field': True,
        'editable': True,
    }
    overrides = {}

    for c, col_data_dict in zip(columns, column_data):
        colmodel = base_colmodel_override.copy()
        colmodel['name'] = c
        colmodel['index'] = c
        colmodel.update(col_data_dict)
        colmodel['search'] = model_field  # Appears in filter dropdown. Currently has to be model due to FK lookup
        colmodel['model_field'] = model_field
        colmodel['queryset_field'] = queryset_field
        overrides[c] = colmodel
    return overrides
