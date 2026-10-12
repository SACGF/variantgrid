from eventlog import views
from eventlog.grids import ApiRequestLogColumns, EventColumns
from library.django_utils import require_superuser
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('', views.eventlog, name='eventlog'),
    path('detail/<int:pk>', views.eventlog_detail, name='eventlog_detail'),
    path('create_event', views.create_event, name='create_event'),
    path('datatable', DatabaseTableView.as_view(column_class=EventColumns), name='event_log_datatable'),
    path('api_request_log', views.api_request_log, name='api_request_log'),
    path('api_request_log/detail/<int:pk>', views.api_request_log_detail, name='api_request_log_detail'),
    path('api_request_log/datatable', require_superuser(DatabaseTableView.as_view(column_class=ApiRequestLogColumns)),
         name='api_request_log_datatable'),
]
