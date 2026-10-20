import logging
from datetime import date
from typing import Any, Optional

from django.contrib.contenttypes.models import ContentType
from django.db.models import QuerySet
from django.http import HttpRequest

from eventlog.models import ApiRequestLog, Event
from library.enums.log_level import LogLevel
from library.guardian_utils import bot_group
from library.utils import emoji_to_unicode
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


class EventColumns(DatatableConfig[Event]):

    def render_data(self, row: dict[str, Any]):
        if filename := row.get('filename'):
            return filename
        elif detail := row.get('details'):
            d = detail.split('\n', 1)[0]
            return emoji_to_unicode(d)

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.set_expand_ajax('eventlog_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('date', client_renderer='TableFormat.timestamp', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('severity', client_renderer='severityRenderer', orderable=True),
            RichColumn('user__username', name='user', orderable=True),
            RichColumn('app_name', label='App name', orderable=True),
            RichColumn('name', orderable=True),
            RichColumn(name='data', extra_columns=["details", "filename"], renderer=self.render_data, client_renderer='TableFormat.limit', client_renderer_kwargs={"limit": 75}),
            RichColumn('id', visible=False)
        ]

    def get_initial_queryset(self) -> QuerySet[Event]:
        event_qs = Event.objects.all()
        if not self.user.is_superuser:
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
                logging.warning("Unexpected filter %s", filter_param)

        exclude_admin = self.get_query_json('exclude_admin')
        if exclude_admin:
            qs = qs.exclude(user__is_superuser=True).exclude(user__isnull=True).exclude(user__groups=bot_group())

        return qs


class ApiRequestLogColumns(DatatableConfig[ApiRequestLog]):
    search_box_enabled = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.set_expand_ajax('api_request_log_detail', expected_height=200)
        self.rich_columns = [
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True, default_sort=SortOrder.DESC),
            self.user_column(),
            RichColumn('method', orderable=True),
            RichColumn('view_name', label='View', orderable=True),
            RichColumn('status_code', label='Status', orderable=True),
            RichColumn('duration_ms', label='ms', orderable=True),
            RichColumn('num_changes', label='Changes', orderable=True, extra_columns=['changes'],
                       renderer=self.render_changes, client_renderer='renderApiRequestLogChanges'),
            RichColumn('path', search=['path', 'request_body'], client_renderer='TableFormat.limit', client_renderer_kwargs={"limit": 75}),
            RichColumn('id', visible=False)
        ]

    @staticmethod
    def render_changes(row: CellData) -> dict[str, Any]:
        changes = row['changes']
        return {
            "recorded": bool(changes),
            "num_changes": row['num_changes'],
            "summary": ApiRequestLog.changes_summary(changes),
        }

    def get_initial_queryset(self) -> QuerySet[ApiRequestLog]:
        if not self.user.is_superuser:
            return ApiRequestLog.objects.none()
        if object_pk := self.get_query_param('object_pk'):
            content_type = ContentType.objects.get_by_natural_key(self.get_query_param('app_label'),
                                                                  self.get_query_param('model_name'))
            return ApiRequestLog.for_object(content_type, object_pk)
        return ApiRequestLog.objects.all()

    def _get_query_date(self, param: str) -> Optional[date]:
        """ An ISO yyyy-mm-dd param (what <input type="date"> submits); anything else is no filter """
        if raw := self.get_query_param(param):
            try:
                return date.fromisoformat(raw)
            except ValueError:
                pass
        return None

    def filter_queryset(self, qs: QuerySet[ApiRequestLog]) -> QuerySet[ApiRequestLog]:
        if view_name := self.get_query_param('view_name'):
            qs = qs.filter(view_name=view_name)
        if self.get_query_json('errors_only'):
            qs = qs.filter(status_code__gte=400)
        if self.get_query_json('wrote_data'):
            qs = qs.filter(num_changes__gt=0)
        # __date is the local (TIME_ZONE) day, so a day picked on the page means that day as the user sees it
        if date_from := self._get_query_date('date_from'):
            qs = qs.filter(created__date__gte=date_from)
        if date_to := self._get_query_date('date_to'):
            qs = qs.filter(created__date__lte=date_to)
        return qs
