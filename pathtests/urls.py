from rest_framework.urlpatterns import format_suffix_patterns

from library.django_utils.jqgrid_view import JQGridView
from pathtests import views, views_autocomplete, views_rest
from pathtests.grids import PathologyTestOrdersGrid, CasesGrid, PathologyTestsGrid
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('clinician_login', views.clinician_login, name='clinician_login'),
    perm_path('clinician_cases/<int:clinician_id>', views.clinician_cases, name='clinician_cases'),
    perm_path('my_cases_tab', views.my_cases_tab, name='my_cases_tab'),
    perm_path('my_cases', views.my_cases, name='my_cases'),
    perm_path('scientist_cases/<int:user_id>', views.scientist_cases, name='scientist_cases'),
    perm_path('follow_scientist/<int:follow_user_id>', views.follow_scientist, name='follow_scientist'),

    # Pathology Tests
    perm_path('', views.pathology_tests, name='pathology_tests'),
    perm_path('cases', views.cases, name='cases'),
    perm_path('pathology_test_requests', views.pathology_test_requests, name='pathology_test_requests'),
    perm_path('manage_pathology_tests', views.manage_pathology_tests, name='manage_pathology_tests'),
    perm_path('view_pathology_test_version/<int:pk>', views.view_pathology_test_version, name='view_pathology_test_version'),
    perm_path('view_pathology_test/<name>', views.view_pathology_test, name='view_pathology_test'),
    perm_path('modify_pathology_test_version/<int:pk>', views.modify_pathology_test_version, name='modify_pathology_test_version'),

    perm_path('view_pathology_test_order/<int:pk>', views.view_pathology_test_order, name='view_pathology_test_order'),
    perm_path('view_external_pathology_test_order/<external_pk>', views.view_external_pathology_test_order, name='view_external_pathology_test_order'),
    perm_path('view_case/<int:pk>', views.view_case, name='view_case'),
    perm_path('view_external_case/<external_pk>', views.view_external_case, name='view_external_case'),
    # Grids
    perm_path('pathology_test_orders/grid/<slug:op>/', JQGridView.as_view(grid=PathologyTestOrdersGrid), name='pathology_test_orders_grid'),
    perm_path('cases/grid/<slug:op>/', JQGridView.as_view(grid=CasesGrid), name='cases_grid'),
    perm_path('pathology_test/grid/<slug:op>/', JQGridView.as_view(grid=PathologyTestsGrid), name='pathology_tests_grid'),
    # Autocompletes
    perm_path('autocomplete/PathologyTest/v2', views_autocomplete.PathologyTestAutocompleteView.as_view(), name='pathology_test_autocomplete'),
    perm_path('autocomplete/PathologyTestVersion', views_autocomplete.PathologyTestVersionAutocompleteView.as_view(), name='pathology_test_version_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/view_pathology_test_version/<int:pk>', views_rest.PathologyTestVersionView.as_view(), name='api_view_pathology_test_version'),
    perm_path('api/view_latest_pathology_test_version/<name>', views_rest.PathologyTestLatestVersionView.as_view(), name='api_view_latest_pathology_test_version'),
]
urlpatterns += format_suffix_patterns(rest_urlpatterns)
