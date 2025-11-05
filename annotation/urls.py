from annotation import views, views_rest
from annotation.grids import VariantAnnotationVersionGrid, AnnotationRunColumns, \
    VariantAnnotationVersionColumns
from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('', views.annotation, name='annotation'),
    path('annotation_build_detail/<genome_build_name>/', views.annotation_build_detail, name='annotation_build_detail'),
    path('annotation_detail', views.annotation_detail, name='annotation_detail'),
    path('annotation_versions', views.annotation_versions, name='annotation_versions'),
    path('version_diffs', views.version_diffs, name='version_diffs'),
    path('view_version_diff/<int:version_diff_id>/', views.view_version_diff, name='view_version_diff'),
    path('variant_annotation_runs', views.variant_annotation_runs, name='variant_annotation_runs'),
    path('annotation_run/view/<int:annotation_run_id>', views.view_annotation_run, name='view_annotation_run'),
    path('annotation_run/retry/<int:annotation_run_id>', views.retry_annotation_run, name='retry_annotation_run'),
    path('annotation_run/retry_upload/<int:annotation_run_id>', views.retry_annotation_run_upload, name='retry_annotation_run_upload'),
    path('create_manual_variant_entry_from_text/<genome_build_name>/<variants_text>', views.create_manual_variant_entry_from_text, name='create_manual_variant_entry_from_text'),
    path('view_annotation_descriptions', views.view_annotation_descriptions, name='view_annotation_descriptions'),
    path('view_annotation_descriptions/<genome_build_name>', views.view_annotation_descriptions, name='view_annotation_descriptions_genome_build'),
    path('about_new_vep_columns', views.about_new_vep_columns, name='about_new_vep_columns'),
    path('view_annotation_version_details/<int:annotation_version_id>', views.view_annotation_version_details, name='view_annotation_version_details'),
    path('load_cached_web_resource/<pk>', views.load_cached_web_resource, name='load_cached_web_resource'),

    path('annotation_version/datatable/<path:genome_build_name>/', DatabaseTableView.as_view(column_class=VariantAnnotationVersionColumns), name='variant_annotation_version_datatable'),
    path('annotation_version/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=VariantAnnotationVersionGrid), name='variant_annotation_version_grid'),

    path('annotation_run/datatables', DatabaseTableView.as_view(column_class=AnnotationRunColumns), name='annotation_run_datatable'),

    path('citations_json/<path:citations_ids_list>', views.citations_json, name='citations_json'),
    path('citation/<str:citation_id>', views.view_citation, name='view_citation'),
    path('citation/<str:citation_id>/detail', views.view_citation_detail, name='view_citation_detail'),

    path('clinvar/<int:clinvar_variation_id>/detail/<int:min_stars>', views.view_clinvar_records_detail, name='view_clinvar_records_detail'),

    path('api/manual_variant_entry_collection/<int:pk>', views_rest.ManualVariantEntryCollectionView.as_view(),
         name='api_manual_variant_entry_collection'),
    path('api/variant_annotation/<genome_build_name>/<variant_string>', views_rest.VariantAnnotationView.as_view(), name='api_variant_annotation')
]
