from mme import views, views_rest
from mme.views import MMEMatchResultColumns
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('api/match', views_rest.MMEMatchView.as_view(), name='mme_api_match'),

    path('classification/<int:classification_id>', views.mme_classification_panel,
         name='mme_classification_panel'),
    path('submission/<int:submission_id>', views.view_mme_submission, name='mme_view_submission'),
    path('submission/<int:submission_id>/submit', views.submit_mme_submission, name='mme_submit'),
    path('match_results/datatable', DatabaseTableView.as_view(column_class=MMEMatchResultColumns),
         name='mme_match_results_datatable'),
]
