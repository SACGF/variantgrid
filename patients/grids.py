from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from django.shortcuts import get_object_or_404

from annotation.models.models_mim_hpo import MIMMorbid, HumanPhenotypeOntology
from annotation.models.models_phenotype_match import PATIENT_OMIM_PATH, \
    PATIENT_HPO_PATH, PATIENT_GENE_SYMBOL_PATH
from genes.models import GeneVersion
from library.django_utils import get_model_fields
from library.jqgrid_abstract_genes_grid import AbstractGenesGrid
from library.jqgrid_user_row_config import JqGridUserRowConfig
from patients.models import PatientRecords, Patient, PatientRecord
from snpdb.models import GenomeBuild


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
            if term_type == 'omim':
                omim = MIMMorbid.objects.get(description=value)
                patient_id_q = Q(**{PATIENT_OMIM_PATH: omim})
            elif term_type == 'hpo':
                hpo = HumanPhenotypeOntology.objects.get(name=value)  # @UndefinedVariable
                patient_id_q = Q(**{PATIENT_HPO_PATH: hpo})
            elif term_type == 'gene':
                # Match to gene_symbol as there may be multiple with that symbol
                patient_id_q = Q(**{PATIENT_GENE_SYMBOL_PATH: value})
            else:
                msg = f"Unknown term type '{term_type}'"
                raise ValueError(msg)

            patient_id_qs = Patient.objects.filter(patient_id_q).values_list("pk", flat=True)
            queryset = queryset.filter(pk__in=patient_id_qs)

        # Add sample_count to queryset
        annotation_kwargs = {"reference_id": StringAgg("specimen__reference_id", ',', distinct=True),
                             "phenotype_matches": StringAgg(PATIENT_HPO_PATH + "__name", '|', distinct=True),
                             "omim": StringAgg(PATIENT_OMIM_PATH + "__description", '|', distinct=True),
                             "genes": StringAgg(PATIENT_GENE_SYMBOL_PATH, '|', distinct=True),
                             "sample_count": Count("sample", distinct=True),
                             "samples": StringAgg("sample__name", ", ", distinct=True)}
        queryset = queryset.annotate(**annotation_kwargs)
        field_names = self.get_field_names() + list(annotation_kwargs.keys())
        self.queryset = queryset.values(*field_names)

        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        EXTRA_COLUMNS = [
            {'index': 'reference_id', 'name': 'reference_id', 'label': 'Specimen ReferenceIDs'},
            {'index': 'phenotype_matches', 'name': 'phenotype_matches', 'label': 'Phenotype Matches', 'classes': 'no-word-wrap', 'formatter': 'phenotypeFormatter'},
            {'index': 'omim', 'name': 'omim', 'label': 'OMIM', 'classes': 'no-word-wrap', 'formatter': 'omimFormatter'},
            {'index': 'genes', 'name': 'genes', 'label': 'Genes', 'classes': 'no-word-wrap', 'formatter': 'geneFormatter'},
            {'index': 'sample_count', 'name': 'sample_count', 'label': '# samples', 'sorttype': 'int', 'width': '30px'},
            {'index': 'samples', 'name': 'samples', 'label': 'Samples'},
        ]
        colmodels.extend(EXTRA_COLUMNS)
        return colmodels


class PatientRecordsGrid(JqGridUserRowConfig):
    model = PatientRecords
    caption = 'PatientRecords'
    fields = ["id", "uploadedpatientrecords__uploaded_file__user__username", "uploadedpatientrecords__uploaded_file__name"]
    colmodel_overrides = {'id': {'width': 40, 'formatter': 'viewPatientRecordsLink'},
                          'uploadedpatientrecords__uploaded_file__user__username': {'label': 'Uploaded by'},
                          'uploadedpatientrecords__uploaded_file__name': {'label': 'Uploaded File Name'}}

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)
        queryset = self.model.objects.filter(uploadedpatientrecords__uploaded_file__user=user)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})


class PatientRecordGrid(JqGridUserRowConfig):
    model = PatientRecord
    caption = 'PatientRecord'
    fields = get_model_fields(PatientRecord)

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


class PatientHPOGenesGrid(AbstractGenesGrid):
    model = GeneVersion
    caption = "Patient HPO genes"

    def __init__(self, user, patient_id, genome_build_name):
        super().__init__(user)

        self.patient = Patient.get_for_user(user, patient_id)
        self.genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.extra_config.update({'sortname': 'name',
                                  'sortorder': 'asc'})

    def get_caption(self):
        return f"{self.caption} ({self.genome_build})"

    def get_column_names(self):
        return ["phenotypes"]

    def get_sql_params_and_columns(self, request):
        sql = """   SELECT DISTINCT genes_geneversion.gene_symbol_id as name,
                    string_agg(DISTINCT annotation_humanphenotypeontology.name, ', ') as phenotypes
                    from annotation_mimgene
                    join genes_gene on (genes_gene.identifier=gene_id)
                    join genes_geneversion on (genes_geneversion.gene_id = genes_gene.identifier)
                    join annotation_phenotypemim on (annotation_phenotypemim.mim_morbid_id=annotation_mimgene.mim_morbid_id)
                    join annotation_humanphenotypeontology on (annotation_humanphenotypeontology.id=annotation_phenotypemim.hpo_id)
                    where annotation_phenotypemim.hpo_id in %s
                    AND genes_geneversion.genome_build_id = %s
                    group by genes_geneversion.gene_symbol_id
        """

        patient_hpo = self.patient.get_hpo_qs()
        if patient_hpo:
            hpo_list = [str(hpo.pk) for hpo in patient_hpo]
        else:
            hpo_list = [-1]  # Invalid
        params = [tuple(hpo_list), self.genome_build.pk]
        return sql, params, None, False

    def get_labels(self):
        return ['Gene ID', "HPOs"]


class PatientMIMGenesGrid(AbstractGenesGrid):
    model = GeneVersion
    caption = "Patient MIM genes"

    def __init__(self, user, patient_id, genome_build_name):
        super().__init__(user)

        self.patient = Patient.get_for_user(user, patient_id)
        self.genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.extra_config.update({'sortname': 'name',
                                  'sortorder': 'asc'})

    def get_caption(self):
        return f"{self.caption} ({self.genome_build})"

    def get_column_names(self):
        return ["mims"]

    def get_sql_params_and_columns(self, request):
        sql = """   SELECT DISTINCT genes_geneversion.gene_symbol_id as name,
                    string_agg(DISTINCT annotation_mimmorbid.description, ', ') as mims
                    from annotation_mimgene
                    join genes_gene on (genes_gene.identifier=gene_id)
                    join genes_geneversion on (genes_geneversion.gene_id = genes_gene.identifier)
                    join annotation_mimmorbid on (annotation_mimmorbid.accession=annotation_mimgene.mim_morbid_id)
                    where annotation_mimmorbid.accession in %s
                    AND genes_geneversion.genome_build_id = %s
                    group by genes_geneversion.gene_symbol_id
        """

        patient_mim = self.patient.get_mim_qs()
        if patient_mim:
            mim_list = [str(mim.pk) for mim in patient_mim]
        else:
            mim_list = [-1]  # Invalid
        params = [tuple(mim_list), self.genome_build.pk]
        return sql, params, None, False

    def get_labels(self):
        return ['Gene ID', "OMIM"]
