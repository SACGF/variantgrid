from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns
from assay_detailed import views
from assay_detailed.views import AssayDetailedRNAColumns
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

router = routers.DefaultRouter()

urlpatterns = [
    path('assays/<int:assay_detailed_rna_pk>/detail', views.view_assay_detailed_rna_detail,
         name='view_assay_detailed_rna_detail')
]

rest_urlpatterns = [
    path('api/assays/datatables/', DatabaseTableView.as_view(column_class=AssayDetailedRNAColumns), name='assay_detailed_rna_datatables')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
