import csv
import json
from typing import Type, Any, Optional, Iterator

from django.core.exceptions import PermissionDenied
from django.http.response import JsonResponse, HttpResponse, StreamingHttpResponse
from django.urls.base import resolve, reverse
from django.utils.text import slugify
from django.views.generic.base import View

from library.utils import nice_class_name, StashFile


class JQGridViewOp:
    CONFIG = "config"
    HANDLER = "handler"
    EDIT = "edit"
    DOWNLOAD = "download"


class JQGridView(View):
    """ Convenience view to quickly create JQGrids. Called in 4 different ways
        depending on the op (see enum above). Add to urls.py like:

        perm_path('analyses/grid/<slug:op>/', JQGridView.as_view(grid=AnalysesGrid, delete=True), name='analyses_grid'),
    """

    grid: Optional[Type[Any]] = None  # JqGridUserRowConfig (or grid initialised w/user, has delete_row method)
    delete_row = False
    csv_download = False  # via request - can also do via JSON

    @staticmethod
    def create_grid_from_request(request, grid_klass, **kwargs):
        kwargs["user"] = request.user
        # extra_filters are commonly used to filter by something custom via JS
        extra_filters = request.GET.get("extra_filters")
        if extra_filters:
            kwargs["extra_filters"] = json.loads(extra_filters)

        return grid_klass(**kwargs)

    @staticmethod
    def export_grid_as_csv(request, grid_klass, basename, **kwargs):

        # TODO: Need to apply default sorting to grid
        grid = JQGridView.create_grid_from_request(request, grid_klass=grid_klass,
                                                   **kwargs)
        return grid_export_request(request, grid, basename)

    def _load_grid(self, request, *args, **kwargs):
        if self.grid is None:
            msg = f"{nice_class_name(self)}.grid not set"
            raise ValueError(msg)
        else:
            grid_klass: Type[Any] = self.grid

        return self.create_grid_from_request(request, grid_klass, **kwargs)

    def get(self, request, *args, **kwargs):
        op = kwargs.pop("op")
        grid_obj = self._load_grid(request, *args, **kwargs)

        if op == JQGridViewOp.CONFIG:
            current_url = resolve(request.path_info).url_name
            if not grid_obj.url:
                kwargs['op'] = JQGridViewOp.HANDLER
                grid_obj.url = reverse(current_url, kwargs=kwargs)
            grid_config = grid_obj.get_config(as_json=False)
            if self.delete_row:
                kwargs['op'] = JQGridViewOp.EDIT
                grid_config["editurl"] = reverse(current_url, kwargs=kwargs)

            return JsonResponse(grid_config)

        if op == JQGridViewOp.HANDLER:
            return HttpResponse(grid_obj.get_json(request), content_type="application/json")

        if op == JQGridViewOp.DOWNLOAD:
            if not self.csv_download:
                msg = f"CSV download requested but 'csv_download' not set on {nice_class_name(self)}"
                raise PermissionDenied(msg)

            try:
                csv_name = grid_obj.csv_name
            except AttributeError:
                try:
                    csv_name = nice_class_name(grid_obj.model)
                except AttributeError:
                    csv_name = nice_class_name(grid_obj)
            return grid_export_request(request, grid_obj, csv_name)

        msg = f"Unknown get op of {op}"
        raise ValueError(msg)

    def post(self, request, *args, **kwargs):
        op = kwargs.pop("op")
        grid_obj = self._load_grid(request, *args, **kwargs)

        if op == JQGridViewOp.EDIT:
            if not self.delete_row:
                msg = f"Edit/delete requested but 'delete_row' not set on {nice_class_name(self)}"
                raise PermissionDenied(msg)

            pk = request.POST['id']
            edit_oper = request.POST["oper"]
            if edit_oper == 'del':
                grid_obj.delete_row(pk)
            else:
                msg = f"Unknown edit operation {edit_oper}"
                raise PermissionDenied(msg)
            return HttpResponse()

        raise ValueError(f"Unknown post op of {op}")


def grid_export_request(request, grid, basename):
    MAX_FILE_NAME_LENGTH = 100  # Shorted as someone got a DDE error opening it in Windows

    request.GET = request.GET.copy()  # Immutable
    request.GET['rows'] = 0  # No pagination
    items = grid.get_items(request)[2]
    colmodels = grid.get_colmodels()
    csv_iterator = grid_export_csv(colmodels, items)
    basename = basename[:MAX_FILE_NAME_LENGTH]
    response = StreamingHttpResponse(csv_iterator, content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{slugify(basename)}.csv"'
    return response


def grid_export_csv(colmodels, items) -> Iterator[str]:
    # If you make the first letters “ID” of a text file
    # Excel incorrectly assumes you are trying to open an SYLK file.
    label_overrides = {"id": "variant_id"}

    pseudo_buffer = StashFile()
    header, labels = colmodel_header_labels(colmodels, label_overrides=label_overrides)
    # Don't use dictwriter as some sample names may be the same
    writer = csv.writer(pseudo_buffer, dialect='excel', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)

    yield pseudo_buffer.value
    for obj in items:
        # labels dict is in same sorted order as header
        row = [obj.get(k) for k in labels]
        writer.writerow(row)
        yield pseudo_buffer.value


def colmodel_header_labels(colmodels, label_overrides=None):
    labels = {}

    header = []
    for c in colmodels:
        name = c['name']
        column_label = c.get("label", name)
        if label_overrides:
            column_label = label_overrides.get(name, column_label)

        labels[name] = column_label
        header.append(column_label)

    return header, labels
