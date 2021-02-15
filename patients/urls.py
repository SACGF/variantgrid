from library.django_utils.jqgrid_view import JQGridView
from patients import views, views_autocomplete
from patients.grids import PatientListGrid, PatientRecordsGrid, PatientRecordGrid, PatientOntologyGenesGrid
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('patient_record_imports', views.patient_record_imports, name='patient_record_imports'),
    perm_path('view_patient_records/<int:patient_records_id>', views.view_patient_records, name='view_patient_records'),
    perm_path('help/import_patient_records_details', views.import_patient_records_details, name='import_patient_records_details'),
    perm_path('example_upload_csv/empty', views.example_upload_csv_empty, name='example_upload_csv_empty'),
    perm_path('example_upload_csv/all', views.example_upload_csv_all, name='example_upload_csv_all'),
    perm_path('example_upload_csv/no_patients', views.example_upload_csv_no_patients, name='example_upload_csv_no_patients'),

    perm_path('patients', views.patients, name='patients'),
    perm_path('patient/create', views.create_patient, name='create_patient'),
    perm_path('view_patient/<int:patient_id>', views.view_patient, name='view_patient'),
    perm_path('view_patient/contact/<int:patient_id>', views.view_patient_contact_tab, name='view_patient_contact_tab'),
    perm_path('view_patient/patient_specimens/<int:patient_id>', views.view_patient_specimens, name='view_patient_specimens'),
    perm_path('view_patient/genes/<int:patient_id>', views.view_patient_genes, name='view_patient_genes'),
    perm_path('view_patient/modifications/<int:patient_id>', views.view_patient_modifications, name='view_patient_modifications'),

    # Attachments
    perm_path('patient_file_upload/<int:patient_id>', views.patient_file_upload, name='patient_file_upload'),
    perm_path('patient_file_delete/<int:pk>', views.patient_file_delete, name='patient_file_delete'),
    perm_path('view_patient_file_attachment/<int:patient_attachment_id>', views.view_patient_file_attachment, name='view_patient_file_attachment'),
    perm_path('view_patient_file_attachment_thumbnail/<int:patient_attachment_id>', views.view_patient_file_attachment_thumbnail, name='view_patient_file_attachment_thumbnail'),

    perm_path('patient_term_matches', views.patient_term_matches, name='patient_term_matches'),
    perm_path('bulk_patient_term/<int:patient_id_offset>', views.bulk_patient_term, name='bulk_patient_term_offset'),
    perm_path('bulk_patient_term', views.bulk_patient_term, name='bulk_patient_term'),
    perm_path('patient_term_approvals/<int:patient_id_offset>', views.patient_term_approvals, name='patient_term_approvals_offset'),
    perm_path('patient_term_approvals', views.patient_term_approvals, name='patient_term_approvals'),
    perm_path('approve_patient_term', views.approve_patient_term, name='approve_patient_term'),
    perm_path('phenotypes_matches', views.phenotypes_matches, name='phenotypes_matches'),

    # Grids
    perm_path('patient/grid/<slug:op>/', JQGridView.as_view(grid=PatientListGrid, delete_row=True), name='patient_grid'),
    perm_path('patient_records/grid/<slug:op>/', JQGridView.as_view(grid=PatientRecordsGrid), name='patient_records_grid'),
    perm_path('patient_record/grid/<int:patient_records_id>/<slug:op>/', JQGridView.as_view(grid=PatientRecordGrid), name='patient_record_grid'),
    perm_path('patient/ontology/genes/grid/<int:patient_id>/<slug:op>/', JQGridView.as_view(grid=PatientOntologyGenesGrid),
              name='patient_ontology_genes_grid'),

    # Autocomplete
    perm_path('autocomplete/Patient/', views_autocomplete.PatientAutocompleteView.as_view(), name='patient_autocomplete'),
    perm_path('autocomplete/Specimen/', views_autocomplete.SpecimenAutocompleteView.as_view(), name='specimen_autocomplete'),
    perm_path('autocomplete/Clinician/', views_autocomplete.ClinicianAutocompleteView.as_view(), name='clinician_autocomplete'),
    perm_path('autocomplete/ExternalPKAutocompleteView', views_autocomplete.ExternalPKAutocompleteView.as_view(), name='external_pk_autocomplete'),
]
