from functools import partial
from typing import Optional

from django.db.models import QuerySet, StringAgg, TextField, Value
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse

from annotation.models.models_phenotype_match import PATIENT_ONTOLOGY_TERM_PATH
from library.utils import JsonDataType
from ontology.grids import AbstractOntologyGenesConfig
from ontology.models import OntologyService, OntologyTerm
from patients.external_references import ExternalReference
from patients.models import Extraction, Patient, PatientRecord, PatientRecords, Specimen
from patients.models_enums import MatchStatus, NucleicAcid, PatientRecordMatchType, TissueStatus
from seqauto.models import SequencingSample
from snpdb.models import Sample
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


def _render_optional_boolean(row: CellData) -> Optional[str]:
    """ affected/consanguineous are nullable - unknown shows as blank rather than 'no' """
    if row.value is None:
        return None
    return "Yes" if row.value else "No"


class PatientListColumns(DatatableConfig[Patient]):
    """ The ontology term columns are aggregated per patient, so filtering to a single term still shows
        every term that patient has """
    search_box_enabled = True
    server_csv_download = True
    # The unfiltered count is over six aggregates, one of them a permission subquery, and only feeds
    # the "(filtered from N total)" text
    count_unfiltered = False

    FILTER_FIELDS = {
        "first_name": "first_name__icontains",
        "last_name": "last_name__icontains",
        "sex": "sex",
    }

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn('id', search=False, visible=False),
            RichColumn('patient_code', label='Patient', orderable=True, extra_columns=['id'],
                       renderer=partial(_render_patient_link, patient_id_column='id'),
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('external_pk__code', label='External ID', orderable=True),
            RichColumn('family_code', label='Family Code', orderable=True),
            RichColumn('phenotype', orderable=True, client_renderer='TableFormat.limit.bind(null, 200)'),
            RichColumn('affected', orderable=True, search=False, renderer=_render_optional_boolean),
            RichColumn('consanguineous', orderable=True, search=False, renderer=_render_optional_boolean),
            RichColumn('reference_id', label='Specimen Reference IDs', orderable=True),
            RichColumn('hpo', label='HPO', orderable=True, search=False, css_class='ontology-terms',
                       client_renderer="ontologyTermsRenderer.bind(null, 'HP')"),
            RichColumn('omim', label='OMIM', orderable=True, search=False, css_class='ontology-terms',
                       client_renderer="ontologyTermsRenderer.bind(null, 'OMIM')"),
            RichColumn('mondo', label='MONDO', orderable=True, search=False, css_class='ontology-terms',
                       client_renderer="ontologyTermsRenderer.bind(null, 'MONDO')"),
            RichColumn('hgnc', label='Genes', orderable=True, search=False, css_class='ontology-terms',
                       client_renderer="ontologyTermsRenderer.bind(null, 'HGNC')"),
            RichColumn('sample_count', label='# Samples', orderable=True, search=False, css_class='num'),
            RichColumn('samples', orderable=True),
            RichColumn('modified', orderable=True, search=False, default_sort=SortOrder.DESC,
                       client_renderer='TableFormat.timestamp'),
            RichColumn('id', name='delete', label='', orderable=False, search=False,
                       renderer=self.render_delete, client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[Patient]:
        """ Samples carry their VCF's permissions rather than their patient's, so the count/names are of the
            samples this user can see rather than of every sample linked to the patient """
        ontology_path = f"{PATIENT_ONTOLOGY_TERM_PATH}__name"

        def ontology_terms(ontology_service: OntologyService) -> StringAgg:
            q_service = Q(**{f"{PATIENT_ONTOLOGY_TERM_PATH}__ontology_service": ontology_service})
            return StringAgg(ontology_path, Value('|'), filter=q_service, distinct=True, output_field=TextField())

        visible_samples_q = Q(sample__in=Sample.filter_for_user(self.user))
        return Patient.filter_for_user(self.user).annotate(
            reference_id=StringAgg("specimen__reference_id", Value(','), distinct=True, output_field=TextField()),
            hpo=ontology_terms(OntologyService.HPO),
            omim=ontology_terms(OntologyService.OMIM),
            mondo=ontology_terms(OntologyService.MONDO),
            hgnc=ontology_terms(OntologyService.HGNC),
            sample_count=Count("sample", filter=visible_samples_q, distinct=True),
            samples=StringAgg("sample__name", Value(", "), filter=visible_samples_q, distinct=True,
                              output_field=TextField()),
        )

    def filter_queryset(self, qs: QuerySet[Patient]) -> QuerySet[Patient]:
        if term_type := self.get_query_param("term_type"):
            term_name = self.get_query_param("term")
            ontology_term = OntologyTerm.objects.get(name=term_name, ontology_service=OntologyService(term_type))
            # Filter to a sub patients list NOT just to certain terms, so that StringAgg will
            # return all the terms for that person
            patient_id_qs = Patient.objects.filter(**{PATIENT_ONTOLOGY_TERM_PATH: ontology_term})
            qs = qs.filter(pk__in=patient_id_qs.values_list("pk", flat=True))

        for param, field in self.FILTER_FIELDS.items():
            if value := self.get_query_param(param):
                qs = qs.filter(**{field: value})
        return qs


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
        self.download_csv_button_enabled = True

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
        self.download_csv_button_enabled = True

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


def _render_extraction_reference(row: CellData) -> JsonDataType:
    """ The claim as the client posted it - a bare reference or the ExternalPK form """
    try:
        return str(ExternalReference.from_data(row.value))
    except ValueError:
        return str(row.value)


def _unmatched_extraction_columns() -> list[RichColumn]:
    """ The parked claim and how long it has been parked - shared by everything that can name an
        extraction before the extraction exists """
    return [
        RichColumn('extraction_reference', label='Claimed reference', orderable=False, search=False,
                   renderer=_render_extraction_reference),
        RichColumn('extraction_match_status', label='Status', orderable=True, search=False,
                   client_renderer=RichColumn.choices_client_renderer(MatchStatus.choices)),
        RichColumn('extraction_match_error', label='Problem', orderable=False),
        RichColumn('extraction_match_date', label='Waiting since', orderable=True, search=False,
                   default_sort=SortOrder.ASC, client_renderer='TableFormat.timestamp'),
    ]


class UnmatchedExtractionSampleColumns(DatatableConfig[Sample]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True
        self.download_csv_button_enabled = True

        self.rich_columns = [
            RichColumn('id', visible=False),
            RichColumn('name', label='Sample', orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn('vcf__name', label='VCF', orderable=True,
                       extra_columns=['vcf_id'], renderer=self.render_vcf_link,
                       client_renderer='TableFormat.linkUrl'),
            *_unmatched_extraction_columns(),
        ]

    @staticmethod
    def render_vcf_link(row: CellData) -> JsonDataType:
        vcf_id = row["vcf_id"]
        return {
            "text": row.value or f"({vcf_id})",
            "url": reverse('view_vcf', kwargs={"vcf_id": vcf_id}),
        }

    def get_initial_queryset(self) -> QuerySet[Sample]:
        return Sample.filter_for_user(self.user).filter(extraction__isnull=True,
                                                        extraction_reference__isnull=False)


class UnmatchedExtractionSequencingSampleColumns(DatatableConfig[SequencingSample]):
    """ The seqauto half - a run's sample sheet can name an extraction before it is accessioned """

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True
        self.download_csv_button_enabled = True

        self.rich_columns = [
            RichColumn('id', visible=False),
            RichColumn('sample_name', label='Sequencing sample', orderable=True),
            RichColumn('sample_sheet__sequencing_run__name', label='Sequencing run', orderable=True,
                       renderer=self.render_sequencing_run_link, client_renderer='TableFormat.linkUrl'),
            *_unmatched_extraction_columns(),
        ]

    @staticmethod
    def render_sequencing_run_link(row: CellData) -> JsonDataType:
        sequencing_run_id = row.value
        return {
            "text": sequencing_run_id,
            "url": reverse('view_sequencing_run', kwargs={"sequencing_run_id": sequencing_run_id}),
        }

    def get_initial_queryset(self) -> QuerySet[SequencingSample]:
        return SequencingSample.objects.filter(extraction__isnull=True,
                                               extraction_reference__isnull=False)


class PatientOntologyGenesConfig(AbstractOntologyGenesConfig):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.patient = Patient.get_for_user(self.user, pk=self.get_query_param("patient_id"))

    def _get_ontology_term_ids(self):
        return self.patient.get_ontology_term_ids()
