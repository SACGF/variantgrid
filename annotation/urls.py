from rest_framework.urlpatterns import format_suffix_patterns

from annotation import views, views_rest
from annotation.grids import VariantAnnotationVersionGrid, AnnotationRunColumns, \
    VariantAnnotationVersionColumns
from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('', views.annotation, name='annotation'),
    perm_path('annotation_build_detail/<genome_build_name>/', views.annotation_build_detail, name='annotation_build_detail'),
    perm_path('annotation_detail', views.annotation_detail, name='annotation_detail'),
    perm_path('annotation_versions', views.annotation_versions, name='annotation_versions'),
    perm_path('version_diffs', views.version_diffs, name='version_diffs'),
    perm_path('view_version_diff/<int:version_diff_id>/', views.view_version_diff, name='view_version_diff'),
    perm_path('variant_annotation_runs', views.variant_annotation_runs, name='variant_annotation_runs'),
    perm_path('annotation_run/view/<int:annotation_run_id>', views.view_annotation_run, name='view_annotation_run'),
    perm_path('annotation_run/retry/<int:annotation_run_id>', views.retry_annotation_run, name='retry_annotation_run'),
    perm_path('annotation_run/retry_upload/<int:annotation_run_id>', views.retry_annotation_run_upload, name='retry_annotation_run_upload'),
    perm_path('create_manual_variant_entry_from_text/<genome_build_name>/<variants_text>', views.create_manual_variant_entry_from_text, name='create_manual_variant_entry_from_text'),
    perm_path('view_annotation_descriptions', views.view_annotation_descriptions, name='view_annotation_descriptions'),
    perm_path('view_annotation_descriptions/<genome_build_name>', views.view_annotation_descriptions, name='view_annotation_descriptions_genome_build'),
    perm_path('about_new_vep_columns', views.about_new_vep_columns, name='about_new_vep_columns'),
    perm_path('view_annotation_version_details/<int:annotation_version_id>', views.view_annotation_version_details, name='view_annotation_version_details'),
    perm_path('load_cached_web_resource/<pk>', views.load_cached_web_resource, name='load_cached_web_resource'),

    perm_path('annotation_version/datatable/<path:genome_build_name>/', DatabaseTableView.as_view(column_class=VariantAnnotationVersionColumns), name='variant_annotation_version_datatable'),
    perm_path('annotation_version/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=VariantAnnotationVersionGrid), name='variant_annotation_version_grid'),

    perm_path('annotation_run/datatables', DatabaseTableView.as_view(column_class=AnnotationRunColumns), name='annotation_run_datatable'),

    perm_path('citations_json/<path:citations_ids_list>', views.citations_json, name='citations_json'),
    perm_path('citation/<str:citation_id>', views.view_citation, name='view_citation'),
    perm_path('citation/<str:citation_id>/detail', views.view_citation_detail, name='view_citation_detail'),

    perm_path('clinvar/<int:clinvar_variation_id>/detail/<int:min_stars>', views.view_clinvar_records_detail, name='view_clinvar_records_detail')
]

rest_urlpatterns = [
    perm_path('api/manual_variant_entry_collection/<int:pk>', views_rest.ManualVariantEntryCollectionView.as_view(),
              name='api_manual_variant_entry_collection'),
    perm_path('api/variant_annotation/<genome_build_name>/<variant_string>', views_rest.VariantAnnotationView.as_view(), name='api_variant_annotation')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
