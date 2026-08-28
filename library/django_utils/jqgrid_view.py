"""
Building a library.jqgrid grid class from a request, and streaming its rows as CSV.

Shared by the DataTables adapter (@see library.django_utils.jqgrid_datatable_adapter) and the
Celery CSV/VCF export.
"""
import json

from library.django_utils.grid_export import csv_streaming_response, grid_export_csv


class JQGridViewOp:
    CONFIG = "config"
    HANDLER = "handler"
    EDIT = "edit"
    DOWNLOAD = "download"


def create_grid_from_request(request, grid_klass, **kwargs):
    kwargs["user"] = request.user
    # extra_filters are commonly used to filter by something custom via JS
    if extra_filters := request.GET.get("extra_filters"):
        kwargs["extra_filters"] = json.loads(extra_filters)

    return grid_klass(**kwargs)


def export_grid_as_csv(request, grid_klass, basename, **kwargs):
    grid = create_grid_from_request(request, grid_klass=grid_klass, **kwargs)
    return grid_export_request(request, grid, basename)


# The variant grids call the column 'id' - say what it's the id of
VARIANT_GRID_LABEL_OVERRIDES = {"id": "variant_id"}


def grid_export_request(request, grid, basename):
    request.GET = request.GET.copy()  # Immutable
    request.GET['rows'] = 0  # No pagination
    items = grid.get_items(request)[2]
    colmodels = grid.get_colmodels()
    csv_iterator = grid_export_csv(colmodels, items, label_overrides=VARIANT_GRID_LABEL_OVERRIDES)
    return csv_streaming_response(basename, csv_iterator)
