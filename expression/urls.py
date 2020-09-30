from expression import views
from expression.grids import ExpressionFilesGrid
from library.django_utils.jqgrid_view import JQGridView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('expression_files', views.expression_files, name='expression_files'),
    perm_path('expression_graph/<int:expression_file_id>', views.expression_graph, name='expression_graph'),
    perm_path('view_expression_file/<int:expression_file_id>', views.view_expression_file, name='view_expression_file'),
    perm_path('expression_files/grid/<slug:op>/', JQGridView.as_view(grid=ExpressionFilesGrid, delete_row=True), name='expression_files_grid'),
]
