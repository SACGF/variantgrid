from django.db.models import QuerySet
from django.http import HttpRequest

from pedigree.models import PedFile, Pedigree
from snpdb.models import UserGridConfig, ImportStatus
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class PedFilesColumns(DatatableConfig[PedFile]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='user__username', label='Uploaded by', orderable=True),
            RichColumn(key='import_status', label='Status', orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[PedFile]:
        return PedFile.filter_for_user(self.user)


class PedigreeColumns(DatatableConfig[Pedigree]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='user__username', label='Created By', orderable=True),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[Pedigree]:
        return Pedigree.filter_for_user(self.user)

    def filter_queryset(self, qs: QuerySet[Pedigree]) -> QuerySet[Pedigree]:
        user_grid_config = UserGridConfig.get(self.user, 'Pedigrees')
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs
