from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns

from manual.views.views import MigrationAttemptColumns, manual_migrations_view
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path

router = routers.DefaultRouter()

urlpatterns = [
    perm_path('manual/migrations', manual_migrations_view, name="manual_migrations"),
    perm_path('manual/migrations/datatable', DatabaseTableView.as_view(column_class=MigrationAttemptColumns), name='manual_migrations_datatable'),
]
rest_urlpatterns = []

urlpatterns += format_suffix_patterns(rest_urlpatterns)
