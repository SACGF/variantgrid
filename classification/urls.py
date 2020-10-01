from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns

from classification.views import views, variant_classification_dashboard_view, \
    variant_classification_export_view, views_autocomplete, variant_classification_import_upload_view
from classification.views.discordance_report_views import discordance_report_view, export_discordance_report
from classification.views.evidence_keys_view import EvidenceKeysView
from classification.views.hgvs_issues_view import view_hgvs_issues, download_hgvs_issues
from classification.views.variant_classification_dashboard_view import problem_download
from classification.views.variant_classification_datatables import VariantClassificationModificationDatatableView
from classification.views.variant_classification_email_view import summary_email_preview_html, \
    summary_email_preview_text
from classification.views.variant_classification_export_view import VariantClassificationApiExportView
from classification.views.variant_classification_overlaps_view import view_overlaps, post_clinical_context, \
    view_clinical_context
from classification.views.variant_classification_view import VariantClassificationView
from classification.views.views import variant_classification_import_tool, AutopopulateView
from variantgrid.perm_path import perm_path


router = routers.DefaultRouter()

urlpatterns = [
    perm_path('activity', views.activity, name='activity'),
    perm_path('activity/<str:latest_timestamp>', views.activity, name='activity_dated'),
    perm_path('classifications', views.classifications, name='classifications'),
    perm_path('create_for_variant/<int:variant_id>', views.create_classification_for_variant, name='create_classification_for_variant'),
    perm_path('create_for_variant_and_transcript/<int:variant_id>/<transcript_id>', views.create_classification_for_variant, name='create_classification_for_variant_and_transcript'),

    # this is uploading the entire import file, distinct from attaching a file to a classification
    perm_path('variant_classification/import_upload', variant_classification_import_upload_view.FileUploadView.as_view(), name="variant_classification_import_upload"),

    perm_path('variant_classification/create', views.create_variant_classification, name='create_variant_classification'),
    perm_path('variant_classification/create_from_hgvs/<genome_build_name>/<hgvs_string>', views.create_classification_from_hgvs, name='create_classification_from_hgvs'),
    perm_path('variant_classification_file_upload/<int:variant_classification_id>', views.variant_classification_file_upload, name='variant_classification_file_upload'),
    perm_path('variant_classification_file_delete/<int:pk>', views.variant_classification_file_delete, name='variant_classification_file_delete'),
    perm_path('view_variant_classification_file_attachment/<int:pk>', views.view_variant_classification_file_attachment, name='view_variant_classification_file_attachment'),
    perm_path('view_variant_classification_file_attachment_thumbnail/<int:pk>', views.view_variant_classification_file_attachment_thumbnail, name='view_variant_classification_file_attachment_thumbnail'),

    perm_path('variant_classification_grid/export/', views.export_variant_classifications_grid, name='export_variant_classifications_grid'),
    perm_path('variant_classification_grid/export_redcap/', views.export_variant_classifications_grid_redcap, name='export_variant_classifications_grid_redcap'),

    perm_path('diff/', views.view_classification_diff, name='variant_classification_diff'),
    perm_path('redcap_data_dictionary.csv', variant_classification_export_view.redcap_data_dictionary, name='redcap_data_dictionary'),
    perm_path('variant_classification/<record_id>/clinvar.xml', variant_classification_export_view.clinvar_xml, name='variant_classification_clinvar_xml'),
    perm_path('variant_classification/<record_id>/variant_classification.csv', variant_classification_export_view.record_csv, name='variant_classification_csv'),
    perm_path('variant_classification/<record_id>/report.html', variant_classification_export_view.template_report, name='view_template_report'),
    perm_path('variant_classification/<record_id>/history', views.variant_classification_history, name='variant_classification_history'),
    perm_path('variant_classification/<record_id>', views.view_variant_classification, name='view_variant_classification'),

    perm_path('evidence_keys/<max_share_level>', views.evidence_keys, name='evidence_keys_max_share_level'),
    perm_path('evidence_keys', views.evidence_keys, name='evidence_keys'),
    perm_path('evidence_keys_logged_in/<max_share_level>', views.evidence_keys_logged_in, name='evidence_keys_logged_in_max_share_level'),
    perm_path('evidence_keys_logged_in', views.evidence_keys_logged_in, name='evidence_keys_logged_in'),

    perm_path('summary_email', summary_email_preview_html, name='summary_email_html'),
    perm_path('summary_email_text', summary_email_preview_text, name='summary_email_text'),
    perm_path('dashboard', variant_classification_dashboard_view.dashboard, name='variant_classification_dashboard'),
    # giving the dashboard a mode has been deprecated, but a lot of refernces still exist
    perm_path('dashboard/all', variant_classification_dashboard_view.dashboard, name="variant_classification_dashboard_all"),
    perm_path('dashboard_download', problem_download, name='variant_classification_dashboard_download'),

    perm_path('variant_classification/discordance_report/<int:report_id>', discordance_report_view, name='discordance_report'),
    perm_path('variant_classification/discordance_report/<int:report_id>/export', export_discordance_report, name='discordance_export'),

    perm_path('export', variant_classification_export_view.export_view, name='variant_classification_export'),
    perm_path('export_redirect', variant_classification_export_view.export_view_redirector, name='variant_classification_export_redirect'),
    perm_path('import', variant_classification_import_tool, name='variant_classification_import_tool'),

    perm_path('clinical_context', post_clinical_context, name='clinical_context'),
    perm_path('overlaps', view_overlaps, name='overlaps'),
    perm_path('clinical_context/<int:pk>', view_clinical_context, name='clinical_context'),
    perm_path('hgvs_issues', view_hgvs_issues, name='hgvs_issues'),
    perm_path('hgvs_issues_download', download_hgvs_issues, name='hgvs_issues_download'),

    perm_path('classification_graphs', views.classification_graphs, name='classification_graphs'),
    perm_path('autocomplete/EvidenceKey/', views_autocomplete.EvidenceKeyAutocompleteView.as_view(), name='evidence_key_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/classifications/auto_populate', AutopopulateView.as_view(), name='variant_classification_auto_populate_api'),

    perm_path('api/classifications/record/', VariantClassificationView.as_view(), name='variant_classification_api'),
    perm_path('api/classifications/record/<record_id>', VariantClassificationView.as_view(), name='variant_classification_with_record_api'),
    perm_path('api/evidence_keys', EvidenceKeysView.as_view(), name='evidence_keys_api'),

    perm_path('api/classifications/v1/record/', VariantClassificationView.as_view(), name='variant_classification_api'),
    perm_path('api/classifications/v1/record/<record_id>', VariantClassificationView.as_view(), name='variant_classification_with_record_api'),

    perm_path('api/classifications/v2/record/', VariantClassificationView.as_view(api_version=2), name='variant_classification_api_2'),
    perm_path('api/classifications/v2/record/<record_id>', VariantClassificationView.as_view(api_version=2), name='variant_classification_with_record_api_2'),

    perm_path('api/classifications/export', VariantClassificationApiExportView.as_view(), name='variant_classification_export_api'),
    perm_path('api/classifications/datatables/', VariantClassificationModificationDatatableView.as_view(), name='variant_classification_datatables')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
