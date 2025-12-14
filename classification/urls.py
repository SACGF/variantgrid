from rest_framework import routers

from classification.views import clinvar_export_view, search_view_metrics
from classification.views import views, classification_dashboard_view, \
    classification_export_view, views_autocomplete, \
    classification_accumulation_graph
from classification.views.allele_grouping_datatables import AlleleGroupingColumns
from classification.views.classification_dashboard_view import issues_download
from classification.views.classification_datatables import ClassificationColumns
from classification.views.classification_email_view import summary_email_preview_html, \
    summary_email_preview_text
from classification.views.classification_export_view import ClassificationApiExportView
from classification.views.classification_grouping_datatables import ClassificationGroupingColumns
from classification.views.classification_overlaps_view import view_overlaps, post_clinical_context, \
    view_clinical_context, view_overlaps_detail
from classification.views.classification_overlaps_vus_view import view_overlaps_vus, view_overlaps_vus_detail
from classification.views.classification_view import ClassificationView, LabGeneClassificationCountsView
from classification.views.classification_view_metrics import view_classification_metrics, \
    view_page_metrics_detail
from classification.views.clinvar_export_view import ClinVarMatchView, clinvar_match_detail
from classification.views.condition_match_test_view import condition_match_test_view, \
    condition_match_test_download_view, condition_obsoletes_view
from classification.views.condition_matching_view import condition_matching_view, condition_matchings_view, \
    ConditionTextColumns, ConditionTextMatchingAPI
from classification.views.discordance_report_triage_view import DiscordanceReportTriageView
from classification.views.classification_my_lab_view import view_my_lab, view_my_lab_detail
from classification.views.discordance_report_views import discordance_report_view, export_discordance_report, \
    discordance_reports_view, discordance_reports_history_detail, discordance_reports_active_detail, \
    discordance_report_review, action_discordance_report_review, discordance_reports_download
from classification.views.evidence_keys_view import EvidenceKeysView
from classification.views.exports.classification_export_formatter_redcap import redcap_data_dictionary
from classification.views.imported_allele_info_view import view_imported_allele_info, ImportedAlleleInfoColumns, \
    view_imported_allele_info_detail, download_allele_info
from classification.views.sample_classification_search_view import sample_classification_search, \
    sample_classification_search_results
from classification.views.views import classification_import_tool, AutopopulateView
from classification.views.views_hgvs_resolution_tool import hgvs_resolution_tool
from classification.views.views_uploaded_classifications_unmapped import UploadedClassificationsUnmappedView, \
    UploadedClassificationsUnmappedColumns, download_classification_unmapped_file, \
    view_uploaded_classification_unmapped, view_uploaded_classification_unmapped_detail, \
    view_uploaded_classification_unmapped_validation_detail, upload_classification_unmapped_download_validation
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

router = routers.DefaultRouter()

urlpatterns = [
    path('activity', views.activity, name='activity'),
    path('activity/lab/<int:lab_id>', views.activity, name='activity_lab'),
    path('activity/user/<int:user_id>', views.activity, name='activity_user'),
    path('activity/report/<int:discordance_report_id>', views.activity, name='activity_discordance'),
    path('classifications', views.classifications, name='classifications'),

    path('groupings', views.classification_groupings, name='classification_groupings'),
    path('groupings/<int:classification_grouping_id>', views.view_classification_grouping_detail, name='classification_grouping_detail'),
    path('groupings/<int:classification_grouping_id>/records', views.view_classification_grouping_records_detail,
         name='classification_grouping_records_detail'),
    path('allele_groupings', views.allele_groupings, name='allele_groupings'),
    path('allele_groupings/<str:lab_id>', views.allele_groupings, name='allele_groupings_lab'),
    path('allele_grouping/<allele_grouping_id>', views.view_allele_grouping_detail, name='allele_grouping_detail'),

    path('create_for_variant/<int:variant_id>/<genome_build_name>', views.CreateClassificationForVariantView.as_view(),
         name='create_classification_for_variant'),

    path('classification/view_metrics', view_classification_metrics, name="classification_view_metrics"),
    path('classification/view_metrics/detail', view_page_metrics_detail, name="classification_view_metrics_detail"),

    # this is uploading the entire import file, distinct from attaching a file to a classification
    path('classification/import_upload', UploadedClassificationsUnmappedView.as_view(), name="classification_upload_unmapped"),
    # TODO move file lab into another subfolder as it gets a bit confused with upload page
    path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>', view_uploaded_classification_unmapped, name="classification_upload_unmapped_status"),
    path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/download_validation', upload_classification_unmapped_download_validation, name="classification_upload_unmapped_status_download_validation"),
    path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/detail', view_uploaded_classification_unmapped_detail, name="classification_upload_unmapped_status_detail"),
    path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/validation_detail', view_uploaded_classification_unmapped_validation_detail, name="classification_upload_unmapped_status_validation_detail"),
    path('classification/import_upload/datatable', DatabaseTableView.as_view(column_class=UploadedClassificationsUnmappedColumns), name='classification_upload_unmapped_datatable'),
    path('classification/import_upload/download/<int:uploaded_classification_unmapped_id>', download_classification_unmapped_file, name='classification_upload_unmapped_download'),
    path('classification/import_upload/<str:lab_id>', UploadedClassificationsUnmappedView.as_view(), name="classification_upload_unmapped_lab"),

    path('classification/create', views.create_classification, name='create_classification'),
    path('classification/create_from_hgvs/<genome_build_name>/<hgvs_string>', views.create_classification_from_hgvs, name='create_classification_from_hgvs'),
    path('classification_file_upload/<int:classification_id>', views.classification_file_upload, name='classification_file_upload'),
    path('classification_file_delete/<int:pk>', views.classification_file_delete, name='classification_file_delete'),
    path('view_classification_file_attachment/<int:pk>', views.view_classification_file_attachment, name='view_classification_file_attachment'),
    path('view_classification_file_attachment_thumbnail/<int:pk>', views.view_classification_file_attachment_thumbnail, name='view_classification_file_attachment_thumbnail'),

    path('classification_grid/export/', views.export_classifications_grid, name='export_classifications_grid'),
    path('classification_grid/export_redcap/', views.export_classifications_grid_redcap, name='export_classifications_grid_redcap'),

    path('clinvar_match', ClinVarMatchView.as_view(), name='clinvar_match'),
    path('clinvar_match/<str:clinvar_key_id>', ClinVarMatchView.as_view(), name='clinvar_match'),
    path('clinvar_match/<str:clinvar_key_id>/match_detail', clinvar_match_detail, name='clinvar_match_detail'),

    path('clinvar_export_summary', clinvar_export_view.clinvar_export_summary, name='clinvar_key_summary'),  # version that lets you pick which clinvarkey if you access to multiple
    path('clinvar_export_summary/<str:clinvar_key_id>', clinvar_export_view.clinvar_export_summary, name='clinvar_key_summary'),
    path('clinvar_export_summary/<str:clinvar_key_id>/refresh', clinvar_export_view.clinvar_export_refresh, name='clinvar_export_refresh'),
    path('clinvar_export_summary/<str:clinvar_key_id>/create_batch', clinvar_export_view.clinvar_export_create_batch, name='clinvar_export_create_batch'),
    path('clinvar_export_summary/<str:clinvar_key_id>/download', clinvar_export_view.clinvar_export_download, name='clinvar_key_summary_export_download'),
    path('clinvar_export/datatable', DatabaseTableView.as_view(column_class=clinvar_export_view.ClinVarExportColumns), name='clinvar_exports_datatables'),
    path('clinvar_export/<int:clinvar_export_id>', clinvar_export_view.clinvar_export_review, name='clinvar_export'),
    path('clinvar_export/<int:clinvar_export_id>/detail', clinvar_export_view.clinvar_export_detail, name='clinvar_export_detail'),
    path('clinvar_export/<int:clinvar_export_id>/history', clinvar_export_view.clinvar_export_history, name='clinvar_export_history'),
    path('clinvar_export_batch/datatable', DatabaseTableView.as_view(column_class=clinvar_export_view.ClinVarExportBatchColumns), name='clinvar_export_batch_datatables'),
    path('clinvar_export_batch/<int:clinvar_export_batch_id>/detail', clinvar_export_view.clinvar_export_batch_detail, name='clinvar_export_batch_detail'),
    path('clinvar_export_batch/<int:clinvar_export_batch_id>/download', clinvar_export_view.clinvar_export_batch_download, name='clinvar_export_batch_download'),

    path('condition_matchings', condition_matchings_view, name='condition_matchings'),
    path('condition_matchings/<str:lab_id>', condition_matchings_view, name='condition_matchings_lab'),
    path('condition_matchings/<str:lab_id>', condition_matchings_view, name='condition_matchings_lab'),
    path('condition_matching/datatable/<str:lab_id>', DatabaseTableView.as_view(column_class=ConditionTextColumns), name='condition_text_datatable'),
    path('condition_matching/<int:pk>', condition_matching_view, name='condition_matching'),

    path('condition_match_test', condition_match_test_view, name='condition_match_test'),
    path('condition_match_test_download', condition_match_test_download_view, name='condition_match_test_download'),

    path('condition_obsoletes', condition_obsoletes_view, name='condition_obsoletes'),

    path('diff/', views.view_classification_diff, name='classification_diff'),
    path('redcap_data_dictionary.csv', redcap_data_dictionary, name='redcap_data_dictionary'),
    path('classification/<classification_id>/classification.csv', classification_export_view.record_csv, name='classification_csv'),
    path('classification/<classification_id>/report.html', classification_export_view.template_report, name='view_template_report'),
    path('classification/<classification_id>/history', views.classification_history, name='classification_history'),
    # classification ID might have a version in it (e.g. a dot)
    path('classification/<classification_id>', views.view_classification, name='view_classification'),

    path('evidence_keys/<max_share_level>', views.evidence_keys, name='evidence_keys_max_share_level'),
    path('evidence_keys', views.evidence_keys, name='evidence_keys'),

    path('sample_classification_search/results', sample_classification_search_results,
         name='sample_classification_search_results'),
    path('sample_classification_search', sample_classification_search, name='sample_classification_search'),

    path('summary_email', summary_email_preview_html, name='summary_email_html'),
    path('summary_email/<str:lab_id>', summary_email_preview_html, name='summary_email_html'),
    path('summary_email_text', summary_email_preview_text, name='summary_email_text'),
    path('summary_email_text/<str:lab_id>', summary_email_preview_text, name='summary_email_html'),

    path('dashboard', classification_dashboard_view.classification_dashboard, name='classification_dashboard'),
    path('dashboard/<str:lab_id>', classification_dashboard_view.classification_dashboard, name='classification_dashboard'),
    path('dashboard_graph/<str:lab_id>', classification_dashboard_view.classification_dashboard_graph_detail, name='classification_dashboard_graph_detail'),
    path('dashboard_classification_special/<str:lab_id>', classification_dashboard_view.classification_dashboard_special_detail, name='classification_dashboard_special_detail'),
    path('dashboard_discordances/<str:lab_id>', classification_dashboard_view.classification_dashboard_classification_discordance_detail, name='classification_dashboard_discordance_detail'),

    path('dashboard_download', issues_download, name='classification_dashboard_download'),
    path('dashboard_download/<str:lab_id>', issues_download, name='classification_dashboard_download'),
    # legacy URL
    path('dashboard/all', classification_dashboard_view.classification_dashboard, name="classification_dashboard_all"),


    path('accumulation_data', classification_accumulation_graph.download_report, name="classification_accumulation_data"),

    path('discordance_reports', discordance_reports_view, name='discordance_reports'),
    path('discordance_reports/<str:lab_id>', discordance_reports_view, name='discordance_reports'),
    path('discordance_reports/<str:lab_id>/history_detail', discordance_reports_history_detail, name='discordance_reports_history_detail'),
    path('discordance_reports/<str:lab_id>/active_detail', discordance_reports_active_detail, name='discordance_reports_active_detail'),
    path('discrodance_reports/<str:lab_id>/download', discordance_reports_download, name='discordance_reports_download'),
    # 'classification' is redundant but there'll be other references to these URLs, so keep the URLs valid
    path('classification/discordance_report/<int:discordance_report_id>', discordance_report_view, name='discordance_report_deprecated'),
    path('classification/discordance_report/<int:discordance_report_id>/review', discordance_report_review, name='discordance_report_review'),
    path('classification/discordance_report/<int:discordance_report_id>/export', export_discordance_report, name='discordance_export_deprecated'),

    path('discordance_report/<int:discordance_report_id>', discordance_report_view, name='discordance_report'),
    path('discordance_report/<int:discordance_report_id>/export', export_discordance_report, name='discordance_export'),
    path('discordance_report_review_action/<int:review_id>', action_discordance_report_review, name='discordance_report_review_action'),

    path('discordance_report_triage/<int:discordance_report_triage_id>', DiscordanceReportTriageView.as_view(), name='discordance_report_triage_detail'),

    path('export_search_data', search_view_metrics.download_search_data, name='export_search_data'),
    path('export', classification_export_view.export_view, name='classification_export'),
    path('internal_lab_download', classification_export_view.internal_lab_download, name='internal_lab_download'),
    path('export_redirect', classification_export_view.export_view_redirector, name='classification_export_redirect'),
    path('import', classification_import_tool, name='classification_import_tool'),

    path('hgvs_resolution_tool', hgvs_resolution_tool, name='hgvs_resolution_tool'),

    path('clinical_context', post_clinical_context, name='clinical_context'),
    path('overlaps', view_overlaps, name='overlaps'),
    path('overlaps/<str:lab_id>', view_overlaps, name='overlaps'),
    path('overlaps_detail/<str:lab_id>', view_overlaps_detail, name='overlaps_detail'),
    path('vus', view_overlaps_vus, name='vus'),
    path('vus/<str:lab_id>', view_overlaps_vus, name='vus'),
    path('vus_detail/<str:lab_id>', view_overlaps_vus_detail, name='vus_detail'),
    path('clinical_context/<int:pk>', view_clinical_context, name='clinical_context'),

    path('imported_allele_info', view_imported_allele_info, name='view_imported_allele_info'),
    path('imported_allele_info/<int:allele_info_id>', view_imported_allele_info_detail, name='view_imported_allele_info_detail'),
    path('imported_allele_info/download', download_allele_info, name='imported_allele_info_download'),

    path('my_lab', view_my_lab, name='my_lab'),
    path('my_lab/<str:lab_id>', view_my_lab, name='my_lab_lab'),
    path('my_lab/detail/<str:lab_id>', view_my_lab_detail, name='my_lab_detail'),

    path('classification_graphs', views.classification_graphs, name='classification_graphs'),
    path('lab_gene_classification_counts', views.lab_gene_classification_counts, name='lab_gene_classification_counts'),
    path('clinical_significance_change_data', views.clin_sig_change_data, name='clinical_significance_change_data'),
    path('autocomplete/EvidenceKey/', views_autocomplete.EvidenceKeyAutocompleteView.as_view(), name='evidence_key_autocomplete'),

    path('api/imported_allele_info/datatables/', DatabaseTableView.as_view(column_class=ImportedAlleleInfoColumns), name='imported_allele_info_datatables'),

    path('api/classifications/auto_populate', AutopopulateView.as_view(), name='classification_auto_populate_api'),

    path('api/classifications/record/', ClassificationView.as_view(api_version=1), name='classification_api'),
    path('api/classifications/record/<record_id>', ClassificationView.as_view(api_version=1), name='classification_with_record_api'),
    path('api/evidence_keys', EvidenceKeysView.as_view(), name='evidence_keys_api'),

    path('api/classifications/v1/record/', ClassificationView.as_view(), name='classification_api'),
    path('api/classifications/v1/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),

    path('api/classifications/v2/record/', ClassificationView.as_view(api_version=2), name='classification_api_2'),
    path('api/classifications/v2/record/<record_id>', ClassificationView.as_view(api_version=2), name='classification_with_record_api_2'),

    path('api/classifications/v3/record/', ClassificationView.as_view(api_version=3), name='classification_api_3'),
    path('api/classifications/v3/record/<record_id>', ClassificationView.as_view(api_version=3),
         name='classification_with_record_api_3'),

    path('api/classifications/export', ClassificationApiExportView.as_view(), name='classification_export_api'),
    path('api/classifications/datatables/', DatabaseTableView.as_view(column_class=ClassificationColumns), name='classification_datatables'),
    path('api/classification/groups/datatables/', DatabaseTableView.as_view(column_class=ClassificationGroupingColumns), name='classification_grouping_datatables'),
    path('api/classification/allele_groups/datatables/<str:lab_id>', DatabaseTableView.as_view(column_class=AlleleGroupingColumns), name='allele_grouping_datatables'),

    path('api/classifications/gene_counts/<lab_id>', LabGeneClassificationCountsView.as_view(),
         name='lab_gene_classification_counts_api'),

    path('api/condition_text_matching/<int:pk>', ConditionTextMatchingAPI.as_view(), name='condition_text_matching_api')
]
