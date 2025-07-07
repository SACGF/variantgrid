from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models import TextField
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.shortcuts import get_object_or_404

from annotation.models.models_phenotype_match import PATIENT_ONTOLOGY_TERM_PATH
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from ontology.grids import AbstractOntologyGenesGrid
from ontology.models import OntologyTerm, OntologyService
from patients.models import PatientRecords, Patient, PatientRecord


class PatientListGrid(JqGridUserRowConfig):
    model = Patient
    caption = 'Patients'
    fields = ["id", "external_pk__code", "family_code", "phenotype", "modified", "affected", "consanguineous"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewPatientLink'}}

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        extra_filters = kwargs.pop("extra_filters", {})
        super().__init__(user)
        queryset = Patient.filter_for_user(user)
        if extra_filters:
            term_type = extra_filters["term_type"]
            value = extra_filters["value"]

            # We need to filter to a sub patients list NOT just to certain terms, so that StringAgg will
            # return all the terms for that person
            if term_type in ('HP', 'OMIM', 'MONDO', 'HGNC'):
                ontology_service = OntologyService(term_type)
                ontology_term = OntologyTerm.objects.get(name=value, ontology_service=ontology_service)
                patient_id_q = Q(**{PATIENT_ONTOLOGY_TERM_PATH: ontology_term})
            else:
                msg = f"Unknown term type '{term_type}'"
                raise ValueError(msg)

            patient_id_qs = Patient.objects.filter(patient_id_q).values_list("pk", flat=True)
            queryset = queryset.filter(pk__in=patient_id_qs)

        ontology_path = f"{PATIENT_ONTOLOGY_TERM_PATH}__name"
        q_hpo = Q(**{f"{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.HPO})
        q_omim = Q(**{f"{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.OMIM})
        q_mondo = Q(**{f"{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.MONDO})
        q_hgnc = Q(**{f"{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": OntologyService.HGNC})
        # Add sample_count to queryset
        annotation_kwargs = {"reference_id": StringAgg("specimen__reference_id", ',',
                                                       distinct=True, output_field=TextField()),
                             "hpo": StringAgg(ontology_path, '|',
                                              filter=q_hpo, distinct=True, output_field=TextField()),
                             "omim": StringAgg(ontology_path, '|',
                                               filter=q_omim, distinct=True, output_field=TextField()),
                             "mondo": StringAgg(ontology_path, '|',
                                                filter=q_mondo, distinct=True, output_field=TextField()),
                             "hgnc": StringAgg(ontology_path, '|',
                                               filter=q_hgnc, distinct=True, output_field=TextField()),
                             "sample_count": Count("sample", distinct=True),
                             "samples": StringAgg("sample__name", ", ", distinct=True, output_field=TextField())}
        queryset = queryset.annotate(**annotation_kwargs)
        field_names = self.get_field_names() + list(annotation_kwargs.keys())
        self.queryset = queryset.values(*field_names)

        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        EXTRA_COLUMNS = [
            {'index': 'reference_id', 'name': 'reference_id', 'label': 'Specimen ReferenceIDs'},
            {'index': 'hpo', 'name': 'hpo', 'label': 'HPO', 'classes': 'no-word-wrap', 'formatter': 'hpoFormatter'},
            {'index': 'omim', 'name': 'omim', 'label': 'OMIM', 'classes': 'no-word-wrap', 'formatter': 'omimFormatter'},
            {'index': 'mondo', 'name': 'mondo', 'label': 'MONDO', 'classes': 'no-word-wrap', 'formatter': 'mondoFormatter'},
            {'index': 'hgnc', 'name': 'hgnc', 'label': 'Genes', 'classes': 'no-word-wrap', 'formatter': 'hgncFormatter'},
            {'index': 'sample_count', 'name': 'sample_count', 'label': '# samples', 'sorttype': 'int', 'width': '30px'},
            {'index': 'samples', 'name': 'samples', 'label': 'Samples'},
        ]
        colmodels.extend(EXTRA_COLUMNS)
        return colmodels


class PatientRecordsGrid(JqGridUserRowConfig):
    model = PatientRecords
    caption = 'PatientRecords'
    fields = ["id", "uploadedpatientrecords__uploaded_file__user__username",
              "uploadedpatientrecords__uploaded_file__name"]
    colmodel_overrides = {'id': {'width': 40, 'formatter': 'viewPatientRecordsLink'},
                          'uploadedpatientrecords__uploaded_file__user__username': {'label': 'Uploaded by'},
                          'uploadedpatientrecords__uploaded_file__name': {'label': 'Uploaded File Name'}}

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)
        queryset = self.model.objects.all()
        if not user.is_superuser:
            queryset = queryset.filter(uploadedpatientrecords__uploaded_file__user=user)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class PatientRecordGrid(JqGridUserRowConfig):
    model = PatientRecord
    caption = 'PatientRecord'
    fields = [
        'id',
        'patient_records__id',
        'record_id',
        'valid',
        'validation_message',
        'sample_id',
        'patient_id',
        'patient__first_name',
        'patient__last_name',
        'patient_match',
        'specimen__reference_id',
        'specimen_match',
        'sample_identifier',
        'sample_name',
        'patient_family_code',
        'patient_first_name',
        'patient_last_name',
        'date_of_birth',
        'date_of_death',
        'sex',
        'affected',
        'consanguineous',
        '_deceased',
        'patient_phenotype',
        'specimen_reference_id',
        'specimen_description',
        'specimen_collected_by',
        'specimen_collection_date',
        'specimen_received_date',
        'specimen_mutation_type',
        'specimen_nucleic_acid_source',
        'specimen_age_at_collection_date'
    ]

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        patient_records_id = kwargs.pop('patient_records_id')
        super().__init__(user)

        patient_records = get_object_or_404(PatientRecords, pk=patient_records_id)
        patient_records.uploaded_file.check_can_view(user)
        queryset = self.model.objects.filter(patient_records=patient_records)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class PatientOntologyGenesGrid(AbstractOntologyGenesGrid):
    def __init__(self, user, patient_id):
        self.patient = Patient.get_for_user(user, pk=patient_id)
        super().__init__()

    def _get_ontology_term_ids(self):
        return self.patient.get_ontology_term_ids()
