from eventlog import views
from eventlog.grids import EventColumns
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path


urlpatterns = [
    path('', views.eventlog, name='eventlog'),
    path('detail/<int:pk>', views.eventlog_detail, name='eventlog_detail'),
    path('create_event', views.create_event, name='create_event'),
    path('datatable', DatabaseTableView.as_view(column_class=EventColumns), name='event_log_datatable')
]