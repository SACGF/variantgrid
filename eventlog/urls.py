from eventlog import views
from eventlog.views import EventLogDatatableView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('', views.eventlog, name='eventlog'),
    perm_path('create_event', views.create_event, name='create_event'),
    perm_path('datatable', EventLogDatatableView.as_view(), name='event_log_datatable')
]
