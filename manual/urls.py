from rest_framework import routers

from library.django_utils import require_superuser
from manual.views.views import MigrationAttemptColumns, manual_migrations_view
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

router = routers.DefaultRouter()

urlpatterns = [
    path('manual/migrations', manual_migrations_view, name="manual_migrations"),
    path('manual/migrations/datatable', require_superuser(DatabaseTableView.as_view(column_class=MigrationAttemptColumns)), name='manual_migrations_datatable'),
]
