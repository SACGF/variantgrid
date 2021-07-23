from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns

from classification.views import views, classification_dashboard_view, \
    classification_export_view, views_autocomplete, classification_import_upload_view, \
    classification_accumulation_graph
from classification.views.condition_match_test_view import condition_match_test_view, condition_match_test_download_view
from classification.views.condition_matching_view import condition_matching_view, condition_matchings_view, \
    ConditionTextColumns, ConditionTextMatchingAPI
from classification.views import clinvar_export_view
from classification.views.discordance_report_views import discordance_report_view, export_discordance_report
from classification.views.evidence_keys_view import EvidenceKeysView
from classification.views.hgvs_issues_view import view_hgvs_issues, download_hgvs_issues, AlleleColumns
from classification.views.classification_dashboard_view import problem_download
from classification.views.classification_datatables import ClassificationDatatableConfig
from classification.views.classification_email_view import summary_email_preview_html, \
    summary_email_preview_text
from classification.views.classification_export_view import ClassificationApiExportView
from classification.views.classification_overlaps_view import view_overlaps, post_clinical_context, \
    view_clinical_context
from classification.views.classification_view import ClassificationView, LabGeneClassificationCountsView
from classification.views.views import classification_import_tool, AutopopulateView
from snpdb.views.datatable_view import DatabasetableView
from variantgrid.perm_path import perm_path


router = routers.DefaultRouter()

urlpatterns = [
    perm_path('activity', views.activity, name='activity'),
    perm_path('activity/<str:latest_timestamp>', views.activity, name='activity_dated'),
    perm_path('classifications', views.classifications, name='classifications'),
    perm_path('create_for_variant/<int:variant_id>/<genome_build_name>', views.CreateClassificationForVariantView.as_view(),
              name='create_classification_for_variant'),

    # this is uploading the entire import file, distinct from attaching a file to a classification
    perm_path('classification/import_upload', classification_import_upload_view.FileUploadView.as_view(), name="classification_import_upload"),

    perm_path('classification/create', views.create_classification, name='create_classification'),
    perm_path('classification/create_from_hgvs/<genome_build_name>/<hgvs_string>', views.create_classification_from_hgvs, name='create_classification_from_hgvs'),
    perm_path('classification_file_upload/<int:classification_id>', views.classification_file_upload, name='classification_file_upload'),
    perm_path('classification_file_delete/<int:pk>', views.classification_file_delete, name='classification_file_delete'),
    perm_path('view_classification_file_attachment/<int:pk>', views.view_classification_file_attachment, name='view_classification_file_attachment'),
    perm_path('view_classification_file_attachment_thumbnail/<int:pk>', views.view_classification_file_attachment_thumbnail, name='view_classification_file_attachment_thumbnail'),

    perm_path('classification_grid/export/', views.export_classifications_grid, name='export_classifications_grid'),
    perm_path('classification_grid/export_redcap/', views.export_classifications_grid_redcap, name='export_classifications_grid_redcap'),

    perm_path('clinvar_export', clinvar_export_view.clinvar_exports_view, name='clinvar_exports'),
    perm_path('clinvar_export/datatable', DatabasetableView.as_view(column_class=clinvar_export_view.ClinVarExportRecordColumns), name='clinvar_exports_datatables'),
    perm_path('clinvar_export/<int:pk>', clinvar_export_view.clinvar_export_review_view, name='clinvar_export'),
    perm_path('clinvar_export/<int:pk>/history', clinvar_export_view.clinvar_export_history_view, name='clinvar_export_history'),
    perm_path('clinvar_export_batch/<int:pk>', clinvar_export_view.clinvar_export_batch_view, name='clinvar_export_batch'),
    perm_path('clinvar_export_batch/<int:pk>/download', clinvar_export_view.clinvar_export_batch_download, name='clinvar_export_batch_download'),

    perm_path('condition_matchings', condition_matchings_view, name='condition_matchings'),
    perm_path('condition_matching/datatable', DatabasetableView.as_view(column_class=ConditionTextColumns), name='condition_text_datatable'),
    perm_path('condition_matching/<int:pk>', condition_matching_view, name='condition_matching'),

    perm_path('condition_match_test', condition_match_test_view, name='condition_match_test'),
    perm_path('condition_match_test_download', condition_match_test_download_view, name='condition_match_test_download'),

    perm_path('diff/', views.view_classification_diff, name='classification_diff'),
    perm_path('redcap_data_dictionary.csv', classification_export_view.redcap_data_dictionary, name='redcap_data_dictionary'),
    perm_path('classification/<record_id>/classification.csv', classification_export_view.record_csv, name='classification_csv'),
    perm_path('classification/<record_id>/report.html', classification_export_view.template_report, name='view_template_report'),
    perm_path('classification/<record_id>/history', views.classification_history, name='classification_history'),
    perm_path('classification/<record_id>', views.view_classification, name='view_classification'),

    perm_path('evidence_keys/<max_share_level>', views.evidence_keys, name='evidence_keys_max_share_level'),
    perm_path('evidence_keys', views.evidence_keys, name='evidence_keys'),
    perm_path('evidence_keys_logged_in/<max_share_level>', views.evidence_keys_logged_in, name='evidence_keys_logged_in_max_share_level'),
    perm_path('evidence_keys_logged_in', views.evidence_keys_logged_in, name='evidence_keys_logged_in'),

    perm_path('summary_email', summary_email_preview_html, name='summary_email_html'),
    perm_path('summary_email_text', summary_email_preview_text, name='summary_email_text'),
    perm_path('dashboard', classification_dashboard_view.dashboard, name='classification_dashboard'),
    # giving the dashboard a mode has been deprecated, but a lot of refernces still exist
    perm_path('dashboard/all', classification_dashboard_view.dashboard, name="classification_dashboard_all"),
    perm_path('dashboard_download', problem_download, name='classification_dashboard_download'),

    perm_path('accumulation_data', classification_accumulation_graph.download_report, name="classification_accumulation_data"),

    perm_path('classification/discordance_report/<int:report_id>', discordance_report_view, name='discordance_report'),
    perm_path('classification/discordance_report/<int:report_id>/export', export_discordance_report, name='discordance_export'),

    perm_path('export', classification_export_view.export_view, name='classification_export'),
    perm_path('export_redirect', classification_export_view.export_view_redirector, name='classification_export_redirect'),
    perm_path('import', classification_import_tool, name='classification_import_tool'),

    perm_path('clinical_context', post_clinical_context, name='clinical_context'),
    perm_path('overlaps', view_overlaps, name='overlaps'),
    perm_path('clinical_context/<int:pk>', view_clinical_context, name='clinical_context'),
    perm_path('hgvs_issues', view_hgvs_issues, name='hgvs_issues'),
    perm_path('hgvs_issues/allele/datatable', DatabasetableView.as_view(column_class=AlleleColumns), name='allele_datatable'),
    perm_path('hgvs_issues_download', download_hgvs_issues, name='hgvs_issues_download'),

    perm_path('classification_graphs', views.classification_graphs, name='classification_graphs'),
    perm_path('lab_gene_classification_counts', views.lab_gene_classification_counts, name='lab_gene_classification_counts'),
    perm_path('clinical_significance_change_data', views.clin_sig_change_data, name='clinical_significance_change_data'),
    perm_path('autocomplete/EvidenceKey/', views_autocomplete.EvidenceKeyAutocompleteView.as_view(), name='evidence_key_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/classifications/auto_populate', AutopopulateView.as_view(), name='classification_auto_populate_api'),

    perm_path('api/classifications/record/', ClassificationView.as_view(), name='classification_api'),
    perm_path('api/classifications/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),
    perm_path('api/evidence_keys', EvidenceKeysView.as_view(), name='evidence_keys_api'),

    perm_path('api/classifications/v1/record/', ClassificationView.as_view(), name='classification_api'),
    perm_path('api/classifications/v1/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),

    perm_path('api/classifications/v2/record/', ClassificationView.as_view(api_version=2), name='classification_api_2'),
    perm_path('api/classifications/v2/record/<record_id>', ClassificationView.as_view(api_version=2), name='classification_with_record_api_2'),

    perm_path('api/classifications/export', ClassificationApiExportView.as_view(), name='classification_export_api'),
    perm_path('api/classifications/datatables/', DatabasetableView.as_view(column_class=ClassificationDatatableConfig), name='classification_datatables'),

    perm_path('api/classifications/gene_counts/<lab_id>', LabGeneClassificationCountsView.as_view(),
              name='lab_gene_classification_counts_api'),

    perm_path('api/condition_text_matching/<int:pk>', ConditionTextMatchingAPI.as_view(), name='condition_text_matching_api')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
