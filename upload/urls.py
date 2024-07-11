from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from upload.grids import UploadPipelineModifiedVariantsGrid, UploadPipelineSkippedAnnotationGrid, \
    UploadStepColumns
from upload.views import views
from upload.views.views import view_upload_step_detail
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('', views.upload, name='upload'),
    perm_path('upload_poll', views.upload_poll, name='upload_poll'),
    perm_path('view_uploaded_file/<int:uploaded_file_id>', views.view_uploaded_file, name='view_uploaded_file'),
    perm_path('view_upload_pipeline/<int:upload_pipeline_id>', views.view_upload_pipeline, name='view_upload_pipeline'),
    perm_path('view_upload_pipeline/warnings_and_errors/<int:upload_pipeline_id>', views.view_upload_pipeline_warnings_and_errors, name='view_upload_pipeline_warnings_and_errors'),
    perm_path('upload_retry_import/<int:upload_pipeline_id>', views.upload_retry_import, name='upload_retry_import'),
    # Grids

    perm_path('upload_pipeline/steps/datatables/', DatabaseTableView.as_view(column_class=UploadStepColumns), name='upload_step_datatables'),
    perm_path('upload_pipeline/step/<int:upload_step_id>', view_upload_step_detail, name='upload_step_detail'),
    perm_path('upload_pipeline/skipped_annotation/grid/<int:upload_pipeline_id>/<slug:op>/',
              JQGridView.as_view(grid=UploadPipelineSkippedAnnotationGrid),
              name='upload_pipeline_skipped_annotation_grid'),
    perm_path('upload_pipeline/modified_variants/grid/<int:upload_pipeline_id>/<slug:op>/', JQGridView.as_view(grid=UploadPipelineModifiedVariantsGrid), name='upload_pipeline_modified_variants_grid'),

    perm_path('view_upload_stats/detail', views.view_upload_stats, name='view_upload_stats_detail'),
    perm_path('vcf_import_info_tags/accept/<int:vcf_import_info_id>', views.accept_vcf_import_info_tag, name='accept_vcf_import_info_tag'),
    perm_path('jfu_upload/', views.jfu_upload, name='jfu_upload'),
    perm_path('jfu_delete/<int:pk>', views.jfu_upload_delete, name='jfu_delete'),
    perm_path('uploaded_file/download/<int:pk>', views.DownloadUploadedFile.as_view(), name='download_uploaded_file'),
]
