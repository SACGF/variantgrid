from eventlog import views
from eventlog.grids import EventColumns
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('', views.eventlog, name='eventlog'),
    perm_path('create_event', views.create_event, name='create_event'),
    perm_path('datatable', DatabaseTableView.as_view(column_class=EventColumns), name='event_log_datatable')
]
