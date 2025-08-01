from typing import Any

from django.db.models import QuerySet
from django.http import HttpRequest

from eventlog.models import Event
from library.enums.log_level import LogLevel
from library.guardian_utils import bot_group
from library.utils import emoji_to_unicode
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class EventColumns(DatatableConfig[Event]):

    def render_data(self, row: dict[str, Any]):
        if filename := row.get('filename'):
            return filename
        elif detail := row.get('details'):
            d = detail.split('\n', 1)[0]
            return emoji_to_unicode(d)

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('eventlog_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('date', client_renderer='TableFormat.timestamp', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('severity', client_renderer='severityRenderer', orderable=True),
            RichColumn('user__username', name='user', orderable=True),
            RichColumn('app_name', label='App name', orderable=True),
            RichColumn('name', orderable=True),
            RichColumn(name='data', extra_columns=["details", "filename"], renderer=self.render_data, client_renderer='TableFormat.limit.bind(null, 75)'),
            RichColumn('id', visible=False)
        ]

    def get_initial_queryset(self) -> QuerySet[Event]:
        event_qs = Event.objects.all()
        if not self.user.is_staff:
            event_qs = event_qs.filter(user=self.user)
        return event_qs

    def filter_queryset(self, qs: QuerySet[Event]) -> QuerySet[Event]:
        filter_param = self.get_query_param('filter')
        if filter_param and filter_param != 'everything':
            if filter_param == 'logins':
                qs = qs.filter(name='login')
            elif filter_param == 'errors':
                qs = qs.filter(severity=LogLevel.ERROR)
            elif filter_param == 'warnings_and_errors':
                qs = qs.filter(severity__in=[LogLevel.ERROR, LogLevel.WARNING])
            elif filter_param == 'events':
                qs = qs.filter(severity=LogLevel.INFO)
            elif filter_param == 'searches':
                qs = qs.filter(name='search')
            else:
                print(f'Unexpected filter {filter_param}')

        exclude_admin = self.get_query_json('exclude_admin')
        if exclude_admin:
            qs = qs.exclude(user__is_superuser=True).exclude(user__isnull=True).exclude(user__groups=bot_group())

        return qs
