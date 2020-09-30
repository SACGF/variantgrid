from guardian.shortcuts import get_objects_for_user

from expression.models import CuffDiffFile
from library.jqgrid_user_row_config import JqGridUserRowConfig


class ExpressionFilesGrid(JqGridUserRowConfig):
    model = CuffDiffFile
    caption = 'Expression Files'
    fields = ["id", "name", "user__username", "import_status", 'sample_1', 'sample_2', 'imported_records']
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "expression-icon",
                                      "url_name": "view_expression_file",
                                      "url_object_column": "id"}},
        'user__username': {'label': 'Uploaded by'}
    }

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)
        queryset = get_objects_for_user(user, 'expression.view_cuffdifffile', accept_global_perms=False)
        self.queryset = queryset.values(*self.get_field_names())
