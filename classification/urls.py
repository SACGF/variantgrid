from abc import ABC
from typing import List, TypeVar, Type, Union, Any, Generic

from django.db import models
from django.http import HttpRequest
from django.urls import register_converter, resolve
from rest_framework import routers
from rest_framework.urlpatterns import format_suffix_patterns
from threadlocals.threadlocals import get_current_request

from classification.views import clinvar_export_view
from classification.views import views, classification_dashboard_view, \
    classification_export_view, views_autocomplete, \
    classification_accumulation_graph
from classification.views.classification_dashboard_view import issues_download
from classification.views.classification_datatables import ClassificationColumns
from classification.views.classification_email_view import summary_email_preview_html, \
    summary_email_preview_text
from classification.views.classification_export_view import ClassificationApiExportView
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
from classification.views.discordance_report_views import discordance_report_view, export_discordance_report, \
    discordance_reports_view, discordance_reports_history_detail, discordance_reports_active_detail, \
    discordance_reports_download
from classification.views.evidence_keys_view import EvidenceKeysView
from classification.views.exports.classification_export_formatter_redcap import redcap_data_dictionary
from classification.views.imported_allele_info_view import view_imported_allele_info, ImportedAlleleInfoColumns, \
    view_imported_allele_info_detail, download_allele_info
from classification.views.views import classification_import_tool, AutopopulateView
from classification.views.views_uploaded_classifications_unmapped import UploadedClassificationsUnmappedView, \
    UploadedClassificationsUnmappedColumns, download_classification_unmapped_file, \
    view_uploaded_classification_unmapped, view_uploaded_classification_unmapped_detail, \
    view_uploaded_classification_unmapped_validation_detail, upload_classification_unmapped_download_validation
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path

router = routers.DefaultRouter()

class PBContext:

    def __init__(self):
        self.patterns = []
        self.segments = []
        self.path = ""

    def push(self, segment: str):
        self.segments.append(segment)
        self.path = "/".join(self.segments)

    def pop(self):
        self.segments.pop()
        self.path = "/".join(self.segments)

    def register_view(self, view, **kwargs):
        self.patterns.append(
            perm_path(self.path, view, **kwargs)
        )


class PB:

    def __init__(self, segment: str):
        self.segment = segment
        self._view = None
        self._name = None
        self._children: List[PB] = []

    def view(self, view=None, name=None) -> 'PB':
        self._view = view
        self._name = name
        return self

    def children(self, *args) -> 'PB':
        self._children += args
        return self

    def _url_patterns(self, context: PBContext):
        context.push(self.segment)
        if _view := self._view:
            context.register_view(_view, name=self._name)
        if self.children:
            for child in self._children:
                child._url_patterns(context)
        context.pop()

    def url_patterns(self):
        context = PBContext()
        self._url_patterns(context=context)
        return context.patterns

"""
Could some information then be passed to the view?
If View is wrapped, could put resolution in view builder, or in view
PB.intParam(param="uploaded_classification_unmapped" resolver=upload_classification_unmapped_id).children(

)
"""

T = TypeVar("T", bound=models.Model)


class SecureModelInstance(Generic[T]):

    def __init__(self,
                 fetcher: Any,
                 param: Any):
        self.fetcher = fetcher
        self.param = param

    def instance(self, request: HttpRequest):
        return self.fetcher.instance(self.param, request)

class SecureModelFetcher(Generic[T]):

    _regex = '[0-9]+'

    def __init__(self,
                 model: Type[T],
                 field: str = "pk",
                 check_read: bool = False,
                 check_lab: Union[str, bool] = False):
        self.model = model
        self.check_read = check_read
        self.check_lab = check_lab

    def __call__(self):
        return self

    @property
    def regex(self):
        return self._regex

    def instance(self, param: Any, request: HttpRequest) -> T:
        qs = self.model.objects.filter(**{self.field:param})
        if request.user.is_superuser:
            return qs.get()
        if self.check_lab:
            from snpdb.models import Lab
            qs = qs.filter(lab__in=Lab.valid_labs_qs(request.user))
        inst = qs.get()
        if self.check_read:
            pass
            # todo check read permission
        return inst




class SecureModelConverter(ABC):
    regex = '[0-9]+'

    def fetcher(self) -> SecureModelFetcher:
        pass

    def to_python(self, value):



urlpatterns = \
    PB("classification").children(
        PB("import_upload").view(UploadedClassificationsUnmappedView.as_view(), "classification_upload_unmapped").children(
            PB("file").children(
               PB("<int:uploaded_classification_unmapped_id>").view(view_uploaded_classification_unmapped, name="classification_upload_unmapped_status").children(
                   PB("download_validation").view(upload_classification_unmapped_download_validation, name="classification_upload_unmapped_status_download_validation"),
                   PB("detail").view(view_uploaded_classification_unmapped_detail, name="classification_upload_unmapped_status_detail"),
                   PB("validation_detail").view(view_uploaded_classification_unmapped_validation_detail, name="classification_upload_unmapped_status_validation_detail")
               )
            ),
            PB("datatable").view(DatabaseTableView.as_view(column_class=UploadedClassificationsUnmappedColumns), name='classification_upload_unmapped_datatable'),
            PB("download").children(
                PB("int:uploaded_classification_unmapped_id").view(download_classification_unmapped_file, name='classification_upload_unmapped_download')
            ),
            PB("str:lab_id").view(UploadedClassificationsUnmappedView.as_view(), name="classification_upload_unmapped_lab")
        )
    ).url_patterns()


urlpatterns += [
    perm_path('activity', views.activity, name='activity'),
    perm_path('activity/lab/<int:lab_id>', views.activity, name='activity_lab'),
    perm_path('activity/user/<int:user_id>', views.activity, name='activity_user'),
    perm_path('activity/report/<int:discordance_report_id>', views.activity, name='activity_discordance'),
    perm_path('classifications', views.classifications, name='classifications'),
    perm_path('create_for_variant/<int:variant_id>/<genome_build_name>', views.CreateClassificationForVariantView.as_view(),
              name='create_classification_for_variant'),

    perm_path('classification/view_metrics', view_classification_metrics, name="classification_view_metrics"),
    perm_path('classification/view_metrics/detail', view_page_metrics_detail, name="classification_view_metrics_detail"),

    # this is uploading the entire import file, distinct from attaching a file to a classification
    # perm_path('classification/import_upload', UploadedClassificationsUnmappedView.as_view(), name="classification_upload_unmapped"),
    # # TODO move file lab into another subfolder as it gets a bit confused with upload page
    # perm_path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>', view_uploaded_classification_unmapped, name="classification_upload_unmapped_status"),
    # perm_path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/download_validation', upload_classification_unmapped_download_validation, name="classification_upload_unmapped_status_download_validation"),
    # perm_path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/detail', view_uploaded_classification_unmapped_detail, name="classification_upload_unmapped_status_detail"),
    # perm_path('classification/import_upload/file/<int:uploaded_classification_unmapped_id>/validation_detail', view_uploaded_classification_unmapped_validation_detail, name="classification_upload_unmapped_status_validation_detail"),
    # perm_path('classification/import_upload/datatable', DatabaseTableView.as_view(column_class=UploadedClassificationsUnmappedColumns), name='classification_upload_unmapped_datatable'),
    # perm_path('classification/import_upload/download/<int:uploaded_classification_unmapped_id>', download_classification_unmapped_file, name='classification_upload_unmapped_download'),
    # perm_path('classification/import_upload/<str:lab_id>', UploadedClassificationsUnmappedView.as_view(), name="classification_upload_unmapped_lab"),

    perm_path('classification/create', views.create_classification, name='create_classification'),
    perm_path('classification/create_from_hgvs/<genome_build_name>/<hgvs_string>', views.create_classification_from_hgvs, name='create_classification_from_hgvs'),
    perm_path('classification_file_upload/<int:classification_id>', views.classification_file_upload, name='classification_file_upload'),
    perm_path('classification_file_delete/<int:pk>', views.classification_file_delete, name='classification_file_delete'),
    perm_path('view_classification_file_attachment/<int:pk>', views.view_classification_file_attachment, name='view_classification_file_attachment'),
    perm_path('view_classification_file_attachment_thumbnail/<int:pk>', views.view_classification_file_attachment_thumbnail, name='view_classification_file_attachment_thumbnail'),

    perm_path('classification_grid/export/', views.export_classifications_grid, name='export_classifications_grid'),
    perm_path('classification_grid/export_redcap/', views.export_classifications_grid_redcap, name='export_classifications_grid_redcap'),

    perm_path('clinvar_match', ClinVarMatchView.as_view(), name='clinvar_match'),
    perm_path('clinvar_match/<str:clinvar_key_id>', ClinVarMatchView.as_view(), name='clinvar_match'),
    perm_path('clinvar_match/<str:clinvar_key_id>/match_detail', clinvar_match_detail, name='clinvar_match_detail'),

    perm_path('clinvar_export_summary', clinvar_export_view.clinvar_export_summary, name='clinvar_key_summary'),  # version that lets you pick which clinvarkey if you access to multiple
    perm_path('clinvar_export_summary/<str:clinvar_key_id>', clinvar_export_view.clinvar_export_summary, name='clinvar_key_summary'),
    perm_path('clinvar_export_summary/<str:clinvar_key_id>/frefresh', clinvar_export_view.clinvar_export_refresh, name='clinvar_export_refresh'),
    perm_path('clinvar_export_summary/<str:clinvar_key_id>/create_batch', clinvar_export_view.clinvar_export_create_batch, name='clinvar_export_create_batch'),
    perm_path('clinvar_export_summary/<str:clinvar_key_id>/download', clinvar_export_view.clinvar_export_download, name='clinvar_key_summary_export_download'),
    perm_path('clinvar_export/datatable', DatabaseTableView.as_view(column_class=clinvar_export_view.ClinVarExportColumns), name='clinvar_exports_datatables'),
    perm_path('clinvar_export/<int:clinvar_export_id>', clinvar_export_view.clinvar_export_review, name='clinvar_export'),
    perm_path('clinvar_export/<int:clinvar_export_id>/detail', clinvar_export_view.clinvar_export_detail, name='clinvar_export_detail'),
    perm_path('clinvar_export/<int:clinvar_export_id>/history', clinvar_export_view.clinvar_export_history, name='clinvar_export_history'),
    perm_path('clinvar_export_batch/datatable', DatabaseTableView.as_view(column_class=clinvar_export_view.ClinVarExportBatchColumns), name='clinvar_export_batch_datatables'),
    perm_path('clinvar_export_batch/<int:clinvar_export_batch_id>/detail', clinvar_export_view.clinvar_export_batch_detail, name='clinvar_export_batch_detail'),
    perm_path('clinvar_export_batch/<int:clinvar_export_batch_id>/download', clinvar_export_view.clinvar_export_batch_download, name='clinvar_export_batch_download'),

    perm_path('condition_matchings', condition_matchings_view, name='condition_matchings'),
    perm_path('condition_matchings/<str:lab_id>', condition_matchings_view, name='condition_matchings_lab'),
    perm_path('condition_matching/datatable/<str:lab_id>', DatabaseTableView.as_view(column_class=ConditionTextColumns), name='condition_text_datatable'),
    perm_path('condition_matching/<int:pk>', condition_matching_view, name='condition_matching'),

    perm_path('condition_match_test', condition_match_test_view, name='condition_match_test'),
    perm_path('condition_match_test_download', condition_match_test_download_view, name='condition_match_test_download'),

    perm_path('condition_obsoletes', condition_obsoletes_view, name='condition_obsoletes'),

    perm_path('diff/', views.view_classification_diff, name='classification_diff'),
    perm_path('redcap_data_dictionary.csv', redcap_data_dictionary, name='redcap_data_dictionary'),
    perm_path('classification/<classification_id>/classification.csv', classification_export_view.record_csv, name='classification_csv'),
    perm_path('classification/<classification_id>/report.html', classification_export_view.template_report, name='view_template_report'),
    perm_path('classification/<classification_id>/history', views.classification_history, name='classification_history'),
    perm_path('classification/<classification_id>', views.view_classification, name='view_classification'),  # classificaiton ID might have a version in it (e.g. a dot)

    perm_path('evidence_keys/<max_share_level>', views.evidence_keys, name='evidence_keys_max_share_level'),
    perm_path('evidence_keys', views.evidence_keys, name='evidence_keys'),

    perm_path('summary_email', summary_email_preview_html, name='summary_email_html'),
    perm_path('summary_email_text', summary_email_preview_text, name='summary_email_text'),


    perm_path('dashboard', classification_dashboard_view.classification_dashboard, name='classification_dashboard'),
    perm_path('dashboard/<str:lab_id>', classification_dashboard_view.classification_dashboard, name='classification_dashboard'),
    perm_path('dashboard_graph/<str:lab_id>', classification_dashboard_view.classification_dashboard_graph_detail, name='classification_dashboard_graph_detail'),
    perm_path('dashboard_download', issues_download, name='classification_dashboard_download'),
    perm_path('dashboard_download/<str:lab_id>', issues_download, name='classification_dashboard_download'),
    # legacy URL
    perm_path('dashboard/all', classification_dashboard_view.classification_dashboard, name="classification_dashboard_all"),


    perm_path('accumulation_data', classification_accumulation_graph.download_report, name="classification_accumulation_data"),

    perm_path('discordance_reports', discordance_reports_view, name='discordance_reports'),
    perm_path('discordance_reports/<str:lab_id>', discordance_reports_view, name='discordance_reports'),
    perm_path('discordance_reports/<str:lab_id>/history_detail', discordance_reports_history_detail, name='discordance_reports_history_detail'),
    perm_path('discordance_reports/<str:lab_id>/active_detail', discordance_reports_active_detail, name='discordance_reports_active_detail'),
    perm_path('discrodance_reports/<str:lab_id>/download', discordance_reports_download, name='discordance_reports_download'),
    # 'classification' is redundant but there'll be other references to these URLs, so keep the URLs valid
    perm_path('classification/discordance_report/<int:discordance_report_id>', discordance_report_view, name='discordance_report_deprecated'),
    perm_path('classification/discordance_report/<int:discordance_report_id>/export', export_discordance_report, name='discordance_export_deprecated'),

    perm_path('discordance_report/<int:discordance_report_id>', discordance_report_view, name='discordance_report'),
    perm_path('discordance_report/<int:discordance_report_id>/export', export_discordance_report, name='discordance_export'),

    perm_path('export', classification_export_view.export_view, name='classification_export'),
    perm_path('export_redirect', classification_export_view.export_view_redirector, name='classification_export_redirect'),
    perm_path('import', classification_import_tool, name='classification_import_tool'),

    perm_path('clinical_context', post_clinical_context, name='clinical_context'),
    perm_path('overlaps', view_overlaps, name='overlaps'),
    perm_path('overlaps/<str:lab_id>', view_overlaps, name='overlaps'),
    perm_path('overlaps_detail/<str:lab_id>', view_overlaps_detail, name='overlaps_detail'),
    perm_path('vus', view_overlaps_vus, name='vus'),
    perm_path('vus/<str:lab_id>', view_overlaps_vus, name='vus'),
    perm_path('vus_detail/<str:lab_id>', view_overlaps_vus_detail, name='vus_detail'),
    perm_path('clinical_context/<int:pk>', view_clinical_context, name='clinical_context'),

    perm_path('imported_allele_info', view_imported_allele_info, name='view_imported_allele_info'),
    perm_path('imported_allele_info/<int:allele_info_id>', view_imported_allele_info_detail, name='view_imported_allele_info_detail'),
    perm_path('imported_allele_info/download', download_allele_info, name='imported_allele_info_download'),

    perm_path('classification_graphs', views.classification_graphs, name='classification_graphs'),
    perm_path('lab_gene_classification_counts', views.lab_gene_classification_counts, name='lab_gene_classification_counts'),
    perm_path('clinical_significance_change_data', views.clin_sig_change_data, name='clinical_significance_change_data'),
    perm_path('autocomplete/EvidenceKey/', views_autocomplete.EvidenceKeyAutocompleteView.as_view(), name='evidence_key_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/imported_allele_info/datatables/', DatabaseTableView.as_view(column_class=ImportedAlleleInfoColumns), name='imported_allele_info_datatables'),

    perm_path('api/classifications/auto_populate', AutopopulateView.as_view(), name='classification_auto_populate_api'),

    perm_path('api/classifications/record/', ClassificationView.as_view(), name='classification_api'),
    perm_path('api/classifications/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),
    perm_path('api/evidence_keys', EvidenceKeysView.as_view(), name='evidence_keys_api'),

    perm_path('api/classifications/v1/record/', ClassificationView.as_view(), name='classification_api'),
    perm_path('api/classifications/v1/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),

    perm_path('api/classifications/v2/record/', ClassificationView.as_view(api_version=2), name='classification_api_2'),
    perm_path('api/classifications/v2/record/<record_id>', ClassificationView.as_view(api_version=2), name='classification_with_record_api_2'),

    perm_path('api/classifications/export', ClassificationApiExportView.as_view(), name='classification_export_api'),
    perm_path('api/classifications/datatables/', DatabaseTableView.as_view(column_class=ClassificationColumns), name='classification_datatables'),

    perm_path('api/classifications/gene_counts/<lab_id>', LabGeneClassificationCountsView.as_view(),
              name='lab_gene_classification_counts_api'),

    perm_path('api/condition_text_matching/<int:pk>', ConditionTextMatchingAPI.as_view(), name='condition_text_matching_api')
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
