from snpdb.views.datatable_view import DatabaseTableView
from upload.grids import (
    UploadPipelineModifiedVariantsColumns,
    UploadPipelineSkippedAnnotationColumns,
    UploadStepColumns,
)
from upload.views import views
from upload.views.views import view_upload_step_detail
from upload.views.views_rest import APIAnnotatedDownloadView, APIFileUploadView, APIUploadStatusView
from variantgrid.perm_path import path

urlpatterns = [
    path('', views.upload, name='upload'),
    path('upload_poll', views.upload_poll, name='upload_poll'),
    path('view_uploaded_file/<int:file_upload_id>', views.view_uploaded_file, name='view_uploaded_file'),
    path('view_upload_pipeline/<int:upload_pipeline_id>', views.view_upload_pipeline, name='view_upload_pipeline'),
    path('view_upload_pipeline/warnings_and_errors/<int:upload_pipeline_id>', views.view_upload_pipeline_warnings_and_errors, name='view_upload_pipeline_warnings_and_errors'),
    path('upload_retry_import/<int:upload_pipeline_id>', views.upload_retry_import, name='upload_retry_import'),
    # Grids

    path('upload_pipeline/steps/datatables/', DatabaseTableView.as_view(column_class=UploadStepColumns), name='upload_step_datatables'),
    path('upload_pipeline/step/<int:upload_step_id>', view_upload_step_detail, name='upload_step_detail'),
    path('upload_pipeline/skipped_annotation/datatable/<int:upload_pipeline_id>/',
         DatabaseTableView.as_view(column_class=UploadPipelineSkippedAnnotationColumns),
         name='upload_pipeline_skipped_annotation_datatable'),
    path('upload_pipeline/modified_variants/datatable/<int:upload_pipeline_id>/',
         DatabaseTableView.as_view(column_class=UploadPipelineModifiedVariantsColumns),
         name='upload_pipeline_modified_variants_datatable'),

    path('view_upload_stats/detail', views.view_upload_stats, name='view_upload_stats_detail'),
    path('vcf_import_info_tags/accept/<int:vcf_import_info_id>', views.accept_vcf_import_info_tag, name='accept_vcf_import_info_tag'),
    path('upload_file/', views.upload_file, name='upload_file'),
    path('upload_file_delete/<int:pk>', views.upload_file_delete, name='upload_file_delete'),
    path('uploaded_file/download/<int:pk>', views.DownloadUploadedFile.as_view(), name='download_uploaded_file'),

    # APIs - Django REST framework
    path('api/v1/file_upload', APIFileUploadView.as_view(), name='api_file_upload'),
    path('api/v1/upload_status/<int:file_upload_id>', APIUploadStatusView.as_view(), name='api_upload_status'),
    path('api/v1/upload_status/sha256/<str:sha256_hash>', APIUploadStatusView.as_view(), name='api_upload_status_sha256'),
    path('api/v1/download/<int:file_upload_id>/<str:export_type>', APIAnnotatedDownloadView.as_view(), name='api_annotated_download'),
    path('api/v1/download/sha256/<str:sha256_hash>/<str:export_type>', APIAnnotatedDownloadView.as_view(), name='api_annotated_download_sha256'),

]
