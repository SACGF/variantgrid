from functools import partial

from django.db.models import QuerySet, StringAgg, TextField, Value
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse

from annotation.models.models_phenotype_match import PATIENT_ONTOLOGY_TERM_PATH
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import JsonDataType
from ontology.grids import AbstractOntologyGenesGrid
from ontology.models import OntologyService, OntologyTerm
from patients.models import Extraction, Patient, PatientRecord, PatientRecords, Specimen
from patients.models_enums import NucleicAcid, PatientRecordMatchType, TissueStatus
from snpdb.models import Sample
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


class PatientListGrid(JqGridUserRowConfig):
    model = Patient
    caption = 'Patients'
    fields = ["id", "patient_code", "external_pk__code", "family_code", "phenotype", "modified",
              "affected", "consanguineous"]
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewPatientLink'},
        'patient_code': {'label': 'Patient Code'},
    }

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
        annotation_kwargs = {"reference_id": StringAgg("specimen__reference_id", Value(','),
                                                       distinct=True, output_field=TextField()),
                             "hpo": StringAgg(ontology_path, Value('|'),
                                              filter=q_hpo, distinct=True, output_field=TextField()),
                             "omim": StringAgg(ontology_path, Value('|'),
                                               filter=q_omim, distinct=True, output_field=TextField()),
                             "mondo": StringAgg(ontology_path, Value('|'),
                                                filter=q_mondo, distinct=True, output_field=TextField()),
                             "hgnc": StringAgg(ontology_path, Value('|'),
                                               filter=q_hgnc, distinct=True, output_field=TextField()),
                             "sample_count": Count("sample", distinct=True),
                             "samples": StringAgg("sample__name", Value(", "), distinct=True, output_field=TextField())}
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


class PatientRecordsColumns(DatatableConfig[PatientRecords]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        # self.expand_client_renderer = DatatableConfig._row_expand_ajax('eventlog_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('id', orderable=True, renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('uploadedpatientrecords__file_upload__created', label="Created",
                       client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('uploadedpatientrecords__file_upload__user__username', label="User", orderable=True),
            RichColumn('uploadedpatientrecords__file_upload__name', orderable=True, label="Filename"),
        ]

    def get_initial_queryset(self) -> QuerySet[PatientRecords]:
        # show_group_data = self.get_query_param("patient_records")
        qs = PatientRecords.objects.all()
        if not self.user.is_superuser:
            qs = qs.filter(uploadedpatientrecords__file_upload__user=self.user)
        return qs


class PatientRecordColumns(DatatableConfig[PatientRecord]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        # self.expand_client_renderer = DatatableConfig._row_expand_ajax('eventlog_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('id', orderable=True, renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('record_id', orderable=True),
            RichColumn('valid', orderable=True),
            RichColumn('validation_message', orderable=True),
            RichColumn('sample_id', orderable=True),
            RichColumn('patient_id', orderable=True),
            RichColumn('patient__patient_code', label='Patient Code', orderable=True),
            RichColumn('patient__first_name', orderable=True),
            RichColumn('patient__last_name', orderable=True),
            RichColumn('patient_match', orderable=True,
                       renderer=partial(self._render_patient_match_type, "patient_match")),
            RichColumn('specimen__reference_id', orderable=True),
            RichColumn('specimen_match', orderable=True,
                       renderer=partial(self._render_patient_match_type, "specimen_match")),
            RichColumn('sample_identifier', orderable=True),
            RichColumn('sample_name', orderable=True),
            RichColumn('patient_family_code', orderable=True),
            RichColumn('patient_code', label='Patient Code (de-identified)', orderable=True),
            RichColumn('patient_first_name', orderable=True),
            RichColumn('patient_last_name', orderable=True),
            RichColumn('date_of_birth', orderable=True, client_renderer='TableFormat.timestamp'),
            RichColumn('date_of_death', orderable=True, client_renderer='TableFormat.timestamp'),
            RichColumn('sex', orderable=True)
        ]

    @staticmethod
    def _render_patient_match_type(column_name, row: CellData):
        match_type = None
        if row.value:
            match_type = PatientRecordMatchType(row.value).label
        return match_type

    def get_initial_queryset(self) -> QuerySet[PatientRecord]:
        patient_records_id = self.get_query_param("patient_records")
        patient_records = get_object_or_404(PatientRecords, pk=patient_records_id)
        patient_records.check_can_view(self.user)
        return PatientRecord.objects.filter(patient_records=patient_records)


def _render_patient_link(row: CellData, patient_id_column: str) -> JsonDataType:
    """ The de-identified patient_code is the text, the same as the patient grid shows """
    patient_id = row[patient_id_column]
    return {
        "text": row.value or f"({patient_id})",
        "url": reverse('view_patient', kwargs={"patient_id": patient_id}),
    }


class SpecimenColumns(DatatableConfig[Specimen]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True

        self.rich_columns = [
            RichColumn('id', visible=False),
            RichColumn('reference_id', label='Reference ID', orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn('patient__patient_code', label='Patient', orderable=True,
                       extra_columns=['patient_id'],
                       renderer=partial(_render_patient_link, patient_id_column='patient_id'),
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('tissue__name', label='Tissue', orderable=True),
            RichColumn('tissue_status', label='Tissue Status', orderable=True, search=False,
                       client_renderer=RichColumn.choices_client_renderer(TissueStatus.choices)),
            RichColumn('collection_date', label='Collected', orderable=True,
                       client_renderer='TableFormat.timestamp'),
            RichColumn('received_date', label='Received', orderable=True,
                       client_renderer='TableFormat.timestamp'),
            RichColumn('external_pk__code', label='External ID', orderable=True),
            RichColumn('extraction_count', label='# Extractions', orderable=True, search=False,
                       css_class='num'),
            RichColumn('modified', orderable=True, search=False, default_sort=SortOrder.DESC,
                       client_renderer='TableFormat.timestamp'),
        ]

    def get_initial_queryset(self) -> QuerySet[Specimen]:
        qs = Specimen.filter_for_user(self.user)
        return qs.annotate(extraction_count=Count("extraction", distinct=True))


class ExtractionColumns(DatatableConfig[Extraction]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True

        self.rich_columns = [
            RichColumn('id', visible=False),
            RichColumn('reference_id', label='Reference ID', orderable=True,
                       renderer=self.render_reference_id, client_renderer='TableFormat.linkUrl'),
            RichColumn('specimen__reference_id', label='Specimen', orderable=True,
                       extra_columns=['specimen_id'],
                       renderer=self.render_specimen_link, client_renderer='TableFormat.linkUrl'),
            RichColumn('specimen__patient__patient_code', label='Patient', orderable=True,
                       extra_columns=['specimen__patient_id'],
                       renderer=partial(_render_patient_link, patient_id_column='specimen__patient_id'),
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('nucleic_acid_source', label='Nucleic Acid', orderable=True, search=False,
                       client_renderer=RichColumn.choices_client_renderer(NucleicAcid.choices)),
            RichColumn('extraction_date', label='Extracted', orderable=True,
                       client_renderer='TableFormat.timestamp'),
            RichColumn('external_pk__code', label='External ID', orderable=True),
            RichColumn('sample_count', label='# Samples', orderable=True, search=False, css_class='num'),
            RichColumn('modified', orderable=True, search=False, default_sort=SortOrder.DESC,
                       client_renderer='TableFormat.timestamp'),
        ]

    @staticmethod
    def render_reference_id(row: CellData) -> JsonDataType:
        """ An extraction's own reference is optional - fall back to the pk, as the preview does """
        extraction_id = row["id"]
        return {
            "text": row.value or f"({extraction_id})",
            "url": reverse('view_extraction', kwargs={"extraction_id": extraction_id}),
        }

    @staticmethod
    def render_specimen_link(row: CellData) -> JsonDataType:
        specimen_id = row["specimen_id"]
        return {
            "text": row.value or f"({specimen_id})",
            "url": reverse('view_specimen', kwargs={"specimen_id": specimen_id}),
        }

    def get_initial_queryset(self) -> QuerySet[Extraction]:
        """ Samples carry their VCF's permissions rather than their patient's, so the count is of the
            samples this user can see rather than of every sample off the extraction """
        qs = Extraction.filter_for_user(self.user)
        visible_samples_q = Q(sample__in=Sample.filter_for_user(self.user))
        return qs.annotate(sample_count=Count("sample", filter=visible_samples_q, distinct=True))


class PatientOntologyGenesGrid(AbstractOntologyGenesGrid):
    def __init__(self, user, patient_id):
        self.patient = Patient.get_for_user(user, pk=patient_id)
        super().__init__()

    def _get_ontology_term_ids(self):
        return self.patient.get_ontology_term_ids()
