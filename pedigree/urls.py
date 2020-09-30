from library.django_utils.jqgrid_view import JQGridView
from pedigree import views, views_autocomplete
from pedigree.grids import PedFilesGrid, PedigreeGrid
from variantgrid.perm_path import perm_path


urlpatterns = [
    perm_path('pedigrees', views.pedigrees, name='pedigrees'),
    perm_path('view_pedigree/<int:pedigree_id>', views.view_pedigree, name='view_pedigree'),
    perm_path('ped_files', views.ped_files, name='ped_files'),
    perm_path('view_ped_file/<int:ped_file_id>', views.view_ped_file, name='view_ped_file'),
    perm_path('pedigree_chart/<int:ped_file_id>', views.pedigree_chart, name='pedigree_chart'),

    perm_path('ped_files/grid/<slug:op>/', JQGridView.as_view(grid=PedFilesGrid, delete_row=True), name='ped_files_grid'),
    perm_path('pedigree/grid/<slug:op>/', JQGridView.as_view(grid=PedigreeGrid, delete_row=True), name='pedigree_grid'),
    perm_path('create_pedigree_from_cohort_and_ped_file_family/<int:cohort_id>/<int:ped_file_family_id>', views.create_pedigree_from_cohort_and_ped_file_family, name='create_pedigree_from_cohort_and_ped_file_family'),

    perm_path('autocomplete/Pedigree/', views_autocomplete.PedigreeAutocompleteView.as_view(), name='pedigree_autocomplete'),
]
