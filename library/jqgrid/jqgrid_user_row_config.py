from library.guardian_utils import check_can_write
from library.jqgrid.jqgrid import JqGrid
from library.utils import nice_class_name
from snpdb.models import UserGridConfig


class JqGridUserRowConfig(JqGrid):
    """ This saves rows per-user in the snpdb.models.UserGridRowConfig table """

    def __init__(self, user):
        super().__init__()
        self.user = user

        if self.caption is None:
            raise ValueError("Need to set 'caption' on JqGridUserRowConfig class - it is used for UserGridRowConfig")

        grid_rows, row_selections = UserGridConfig.get_rows_and_selections(user, self.caption)
        self.grid_rows = grid_rows
        self.extra_config["autoencode"] = True
        self.extra_config.update({'rowNum': grid_rows, 'rowList': row_selections})

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        for cm in colmodels:
            hide_non_admin = cm.get("hide_non_admin", False)
            if hide_non_admin and not self.user.is_superuser:
                cm["hidden"] = True
        return colmodels

    @property
    def csv_name(self):
        return nice_class_name(self.model)

    def delete_row(self, pk):
        """ Deletes a row record - called by JQGridView if allow_delete is set
            This must be overridden if object does not provide a can_write method """

        obj = self.model.objects.get(pk=pk)
        try:
            check_can_write(obj, self.user)
        except AttributeError as exc:
            msg = "delete_row should be overridden where object does not provide a can_write method "
            raise NotImplementedError(msg) from exc
        obj.delete()
