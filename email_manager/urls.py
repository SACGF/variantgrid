from email_manager.views.email_manager_views import email_manager_view, EmailColumns, email_detail, email_pure
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('', email_manager_view, name='email_log'),
    path('detail/<int:email_id>', email_detail, name='email_detail'),
    path('pure/<int:email_id>', email_pure, name='email_pure'),
    path('datatable', DatabaseTableView.as_view(column_class=EmailColumns), name='email_datatable')
]
