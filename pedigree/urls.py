from pedigree import views, views_autocomplete
from pedigree.grids import PedFilesColumns, PedigreeColumns
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path


urlpatterns = [
    path('pedigrees', views.pedigrees, name='pedigrees'),
    path('view_pedigree/<int:pedigree_id>', views.view_pedigree, name='view_pedigree'),
    path('ped_files', views.ped_files, name='ped_files'),
    path('view_ped_file/<int:ped_file_id>', views.view_ped_file, name='view_ped_file'),
    path('pedigree_chart/<int:ped_file_id>', views.pedigree_chart, name='pedigree_chart'),

    path('ped_files/datatables/', DatabaseTableView.as_view(column_class=PedFilesColumns), name='ped_files_datatables'),
    path('pedigree/datatables/', DatabaseTableView.as_view(column_class=PedigreeColumns), name='pedigree_datatables'),
    path('create_pedigree_from_cohort_and_ped_file_family/<int:cohort_id>/<int:ped_file_family_id>', views.create_pedigree_from_cohort_and_ped_file_family, name='create_pedigree_from_cohort_and_ped_file_family'),

    path('autocomplete/Pedigree/', views_autocomplete.PedigreeAutocompleteView.as_view(), name='pedigree_autocomplete'),
]
