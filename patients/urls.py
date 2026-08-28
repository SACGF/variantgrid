from django.urls import include
from rest_framework import routers

from library.django_utils.datatable_dataframe import DataFrameTableView
from patients import views, views_autocomplete
from patients.grids import (
    ExtractionColumns,
    PatientListColumns,
    PatientOntologyGenesConfig,
    PatientRecordColumns,
    PatientRecordsColumns,
    SpecimenColumns,
    UnmatchedExtractionSampleColumns,
    UnmatchedExtractionSequencingSampleColumns,
)
from patients.views_rest import (
    ExtractionViewSet,
    PatientViewSet,
    SpecimenMeasureBulkCreateView,
    SpecimenMeasureViewSet,
    SpecimenViewSet,
)
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('patient_imports', views.patient_imports, name='patient_imports'),
    path('patient_import/<int:patient_records_id>', views.view_patient_import, name='view_patient_import'),
    path('patient_import/view_patient_record/<int:pk>', views.view_patient_record, name='view_patient_record'),
    path('help/import_patient_records_details', views.import_patient_records_details, name='import_patient_records_details'),
    path('example_upload_csv/empty', views.example_upload_csv_empty, name='example_upload_csv_empty'),
    path('example_upload_csv/all', views.example_upload_csv_all, name='example_upload_csv_all'),
    path('example_upload_csv/no_patients', views.example_upload_csv_no_patients, name='example_upload_csv_no_patients'),

    path('patients', views.patients, name='patients'),
    path('patient/create', views.create_patient, name='create_patient'),
    path('view_patient/<int:patient_id>', views.view_patient, name='view_patient'),
    path('view_patient/contact/<int:patient_id>', views.view_patient_contact_tab, name='view_patient_contact_tab'),
    path('view_patient/patient_specimens/<int:patient_id>', views.view_patient_specimens, name='view_patient_specimens'),
    path('view_patient/patient_extractions/<int:patient_id>', views.view_patient_extractions, name='view_patient_extractions'),
    path('specimens', views.specimens, name='specimens'),
    path('extractions', views.extractions, name='extractions'),
    path('view_specimen/<int:specimen_id>', views.view_specimen, name='view_specimen'),
    path('view_extraction/<int:extraction_id>', views.view_extraction, name='view_extraction'),
    path('extraction/<int:extraction_id>/samples', views.extraction_samples, name='extraction_samples'),
    path('unmatched_extractions', views.unmatched_extractions, name='unmatched_extractions'),
    path('view_patient/genes/<int:patient_id>', views.view_patient_genes, name='view_patient_genes'),
    path('view_patient/modifications/<int:patient_id>', views.view_patient_modifications, name='view_patient_modifications'),

    # Attachments
    path('patient_file_upload/<int:patient_id>', views.patient_file_upload, name='patient_file_upload'),
    path('patient_file_delete/<int:pk>', views.patient_file_delete, name='patient_file_delete'),
    path('view_patient_file_attachment/<int:patient_attachment_id>', views.view_patient_file_attachment, name='view_patient_file_attachment'),
    path('view_patient_file_attachment_thumbnail/<int:patient_attachment_id>', views.view_patient_file_attachment_thumbnail, name='view_patient_file_attachment_thumbnail'),

    path('patient_term_matches', views.patient_term_matches, name='patient_term_matches'),
    path('bulk_patient_term/<int:patient_id_offset>', views.bulk_patient_term, name='bulk_patient_term_offset'),
    path('bulk_patient_term', views.bulk_patient_term, name='bulk_patient_term'),
    path('patient_term_approvals/<int:patient_id_offset>', views.patient_term_approvals, name='patient_term_approvals_offset'),
    path('patient_term_approvals', views.patient_term_approvals, name='patient_term_approvals'),
    path('approve_patient_term', views.approve_patient_term, name='approve_patient_term'),
    path('phenotypes_matches', views.phenotypes_matches, name='phenotypes_matches'),

    # Grids
    path('patient/datatables/', DatabaseTableView.as_view(column_class=PatientListColumns),
         name='patient_datatables'),
    path('patient_records/datatables/', DatabaseTableView.as_view(column_class=PatientRecordsColumns),
         name='patient_records_datatables'),
    path('patient_record/datatables/', DatabaseTableView.as_view(column_class=PatientRecordColumns),
         name='patient_record_datatables'),
    path('specimen/datatables/', DatabaseTableView.as_view(column_class=SpecimenColumns),
         name='specimen_datatables'),
    path('extraction/datatables/', DatabaseTableView.as_view(column_class=ExtractionColumns),
         name='extraction_datatables'),
    path('unmatched_extraction/sample/datatables/',
         DatabaseTableView.as_view(column_class=UnmatchedExtractionSampleColumns),
         name='unmatched_extraction_sample_datatables'),
    path('unmatched_extraction/sequencing_sample/datatables/',
         DatabaseTableView.as_view(column_class=UnmatchedExtractionSequencingSampleColumns),
         name='unmatched_extraction_sequencing_sample_datatables'),
    path('patient/ontology/genes/datatable/<int:patient_id>/', DataFrameTableView.as_view(column_class=PatientOntologyGenesConfig),
         name='patient_ontology_genes_datatable'),

    # Autocomplete
    path('autocomplete/Patient/', views_autocomplete.PatientAutocompleteView.as_view(), name='patient_autocomplete'),
    path('autocomplete/Specimen/', views_autocomplete.SpecimenAutocompleteView.as_view(), name='specimen_autocomplete'),
    path('autocomplete/Extraction/', views_autocomplete.ExtractionAutocompleteView.as_view(), name='extraction_autocomplete'),
    path('autocomplete/Clinician/', views_autocomplete.ClinicianAutocompleteView.as_view(), name='clinician_autocomplete'),
    path('autocomplete/ExternalPKAutocompleteView', views_autocomplete.ExternalPKAutocompleteView.as_view(), name='external_pk_autocomplete'),
]

router = routers.DefaultRouter()
router.register(r'api/v1/patient', PatientViewSet, basename='api_patient')
router.register(r'api/v1/specimen', SpecimenViewSet, basename='api_specimen')
router.register(r'api/v1/extraction', ExtractionViewSet, basename='api_extraction')
router.register(r'api/v1/specimen_measure', SpecimenMeasureViewSet, basename='api_specimen_measure')

urlpatterns += [
    # Ahead of the router, whose detail route would otherwise read 'bulk_create' as a primary key
    path('api/v1/specimen_measure/bulk_create', SpecimenMeasureBulkCreateView.as_view(),
         name='api_specimen_measure_bulk_create'),
    path('', include(router.urls), name='patients_apis'),
]
