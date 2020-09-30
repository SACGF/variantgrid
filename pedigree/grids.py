from library.jqgrid_user_row_config import JqGridUserRowConfig
from pedigree.models import PedFile, Pedigree


class PedFilesGrid(JqGridUserRowConfig):
    model = PedFile
    caption = '.ped Files'
    fields = ["id", "name", "user__username", 'import_status']
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "pedigree-icon",
                                      "url_name": "view_ped_file",
                                      "url_object_column": "id"}},
        'user__username': {'label': 'Uploaded by'},
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.filter_for_user(user)
        self.queryset = queryset.values(*self.get_field_names())


class PedigreeGrid(JqGridUserRowConfig):
    model = Pedigree
    caption = 'Pedigrees'
    fields = ["id", "name", "user__username"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "pedigree-icon",
                                      "url_name": "view_pedigree",
                                      "url_object_column": "id"}},
        'user__username': {'label': 'Created By'},
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.filter_for_user(user)
        self.queryset = queryset.values(*self.get_field_names())
