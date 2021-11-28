import json

from django.core.exceptions import PermissionDenied
from django.http.response import JsonResponse, HttpResponse
from django.urls.base import resolve, reverse
from django.views.generic.base import View

from library.jqgrid_export import grid_export_request
from library.utils import nice_class_name


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

    grid = None  # JqGridUserRowConfig (or grid initialised w/user, has delete_row method)
    delete_row = False
    csv_download = False  # via request - can also do via JSON

    def _load_grid(self, request, *args, **kwargs):
        grid_klass = self.grid
        if grid_klass is None:
            msg = f"{nice_class_name(self)}.grid not set"
            raise ValueError(msg)

        kwargs["user"] = request.user
        # extra_filters are commonly used to filter by something custom via JS
        extra_filters = request.GET.get("extra_filters")
        if extra_filters:
            kwargs["extra_filters"] = json.loads(extra_filters)

        return grid_klass(**kwargs)

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
                msg = "Edit/delete requested but 'delete_row' not set on %s" % nice_class_name(self)
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
