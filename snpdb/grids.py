import operator
from functools import reduce
from typing import Optional, Any

from django.conf import settings
from django.contrib.postgres.aggregates.general import StringAgg
from django.db.models import F, QuerySet, Value, Func
from django.db.models.aggregates import Count, Max
from django.db.models.fields import CharField, TextField
from django.db.models.query_utils import Q
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse
from guardian.shortcuts import get_objects_for_user

from annotation.models import ManualVariantEntryCollection, PATIENT_ONTOLOGY_TERM_PATH
from genes.models import GeneSymbol, SampleGeneList
from library.django_utils import get_url_from_view_path
from library.genomics.vcf_enums import INFO_LIFTOVER_SWAPPED_REF_ALT
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.unit_percent import get_allele_frequency_formatter
from library.utils import calculate_age, JsonDataType
from ontology.models import OntologyService, OntologyTerm
from snpdb.grid_columns.custom_columns import get_variantgrid_extra_annotate
from snpdb.models import VCF, Cohort, Sample, ImportStatus, \
    GenomicIntervalsCollection, CustomColumnsCollection, Variant, Trio, UserGridConfig, GenomeBuild, ClinGenAllele, \
    VariantZygosityCountCollection, TagColorsCollection, LiftoverRun, AlleleConversionTool, AlleleLiftover, \
    ProcessingStatus, Allele
from snpdb.sample_filters import get_sample_ontology_q, get_sample_qc_gene_list_gene_symbol_q
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs, remove_soft_deleted_vcfs_task
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from uicore.templatetags.js_tags import jsonify_for_js


class VCFListGrid(JqGridUserRowConfig):
    model = VCF
    caption = 'VCFs'
    fields = ["id", "name", "vcf_url", "date", "import_status", "genome_build__name", "user__username", "source",
              "uploadedvcf__uploaded_file__import_source", "genotype_samples", "project__name", "cohort__import_status",
              "uploadedvcf__vcf_importer__name", 'uploadedvcf__vcf_importer__version']
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'width': 550,
                 'formatter': 'viewVCFLink',
                 'formatter_kwargs': {"icon_css_class": "vcf-icon",
                                      "url_name": "view_vcf",
                                      "url_object_column": "id"}},
        "vcf_url": {'name': 'vcf_url', 'label': 'VCF URL', "model_field": False, 'hidden': True},
        'import_status': {'formatter': 'viewImportStatus'},
        "genome_build__name": {"label": "Genome Build"},
        'user__username': {'label': 'Uploaded by', 'width': 60},
        'source': {'label': 'VCF source'},
        "project__name": {'label': "Project"},
        'cohort__import_status': {'hidden': True},
        'uploadedvcf__vcf_importer__name': {"label": 'VCF Importer', "hide_non_admin": True},
        'uploadedvcf__vcf_importer__version': {"label": 'VCF Importer Version', "hide_non_admin": True},
    }

    def __init__(self, user, **kwargs):
        extra_filters = kwargs.get("extra_filters")
        super().__init__(user)
        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = VCF.filter_for_user(user, group_data=user_grid_config.show_group_data)

        # Set via vcf_grid_filter_tags
        if extra_filters:
            if project := extra_filters.get("project"):
                queryset = queryset.filter(project=project)
            if genome_build_name := extra_filters.get("genome_build_name"):
                genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
                queryset = queryset.filter(genome_build=genome_build)

        fake_number = "1234567890"
        view_vcf_url = reverse('view_vcf', kwargs={"vcf_id": fake_number}).rstrip(fake_number)
        view_vcf_url_prefix = get_url_from_view_path(view_vcf_url)
        annotation_kwargs = {
            "vcf_url": Func(
                Value(view_vcf_url_prefix),
                F("pk"),
                function="CONCAT",
                output_field=CharField(),
            ),
        }
        queryset = queryset.annotate(**annotation_kwargs)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, pk):
        """ Do async as it may be slow """
        soft_delete_vcfs(self.user, pk)


# TODO: Merge this an cohort grid below into 1
class SamplesListGrid(JqGridUserRowConfig):
    model = Sample
    caption = 'Samples'
    fields = ["id", "name", "sample_url", "het_hom_count", "vcf__date", "import_status",
              "vcf__genome_build__name", "variants_type", "vcf__user__username", "vcf__source", "vcf__name", "vcf_url",
              "vcf__project__name", "vcf__uploadedvcf__uploaded_file__import_source",
              "sample_gene_list_count", "activesamplegenelist__id",
              "mutationalsignature__id", "mutationalsignature__summary",
              "somaliersampleextract__somalierancestry__predicted_ancestry",
              "patient__first_name", "patient__last_name", "patient__sex",
              "patient__date_of_birth", "patient__date_of_death",
              "specimen__reference_id", "specimen__tissue__name", "specimen__collection_date", "vcf__id"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {"width": 400,
                 'formatter': 'viewSampleLink',
                 'formatter_kwargs': {"icon_css_class": "sample-icon",
                                      "url_name": "view_sample",
                                      "url_object_column": "id"}},
        'import_status': {'formatter': 'viewImportStatus'},
        'vcf__id': {"hidden": True},
        "vcf__genome_build__name": {"label": "Genome Build"},
        'vcf__source': {'label': 'VCF source'},
        'vcf__name': {
            'label': 'VCF Name', "width": 600,
            "formatter": 'linkFormatter',
            'formatter_kwargs': {"icon_css_class": "vcf-icon",
                                 "url_name": "view_vcf",
                                 "url_object_column": "vcf__id"}
        },
        "vcf__project__name": {'label': "Project"},
        "sample_gene_list_count": {'name': 'sample_gene_list_count', 'label': '# Sample GeneLists',
                                   "model_field": False, "formatter": "viewSampleGeneList", 'sorttype': 'int'},
        'activesamplegenelist__id': {'hidden': True},
        'mutationalsignature__id': {'hidden': True},
        'mutationalsignature__summary': {'label': 'Mutational Signature',
                                         'formatter': 'viewMutationalSignature'},
        "somaliersampleextract__somalierancestry__predicted_ancestry": {"label": "Predicted Ancestry"},
        'patient__last_name': {'label': 'Last Name'},
        'patient__sex': {'label': 'Sex'},
        'patient__date_of_birth': {'label': 'D.O.B.'},
        'patient__date_of_death': {'hidden': True},
        'het_hom_count': {'name': 'het_hom_count', "model_field": False, 'sorttype': 'int',
                          'label': 'Het/Hom Count'},
        "specimen__reference_id": {'label': 'Specimen'},
        # These urls are only there for CSV export
        "sample_url": {'name': 'sample_url', 'label': 'Sample URL', "model_field": False, 'hidden': True},
        "vcf_url": {'name': 'vcf_url', 'label': 'VCF URL', "model_field": False, 'hidden': True},
    }

    def __init__(self, user, **kwargs):
        extra_filters = kwargs.get("extra_filters")
        super().__init__(user)

        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = Sample.filter_for_user(user, group_data=user_grid_config.show_group_data)

        # Set via vcf_grid_filter_tags
        if extra_filters:
            if project := extra_filters.get("project"):
                queryset = queryset.filter(vcf__project=project)
            if genome_build_name := extra_filters.get("genome_build_name"):
                genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
                queryset = queryset.filter(vcf__genome_build=genome_build)
            variants_type = extra_filters.get("variants_type")
            if variants_type is not None:
                queryset = queryset.filter(variants_type__in=variants_type)

        # If you don't have permission to view a patient - blank it out
        # If you have read only and
        # TODO: We need to pass whole row in - as we need date of death to display age
        if settings.PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB:
            dob_colmodel = self._overrides.get('patient__date_of_birth', {})
            dob_colmodel['label'] = "Age"
            dob_colmodel['server_side_formatter'] = lambda row, field: calculate_age(row[field])
            self._overrides['patient__date_of_birth'] = dob_colmodel

        # Only show mut sig column if we have any
        if not queryset.filter(mutationalsignature__isnull=False).exists():
            mut_sig_colmodel = self._overrides.get('mutationalsignature__summary', {})
            mut_sig_colmodel['hidden'] = True
            self._overrides['mutationalsignature__summary'] = mut_sig_colmodel

        if not queryset.filter(somaliersampleextract__somalierancestry__isnull=False).exists():
            somalier_ancestry_colmodel = self._overrides.get('somaliersampleextract__somalierancestry__predicted_ancestry', {})
            somalier_ancestry_colmodel['hidden'] = True
            self._overrides['somaliersampleextract__somalierancestry__predicted_ancestry'] = somalier_ancestry_colmodel

        if not queryset.filter(samplegenelist__isnull=False).exists():
            sample_gene_list_count = self._overrides.get('sample_gene_list_count', {})
            sample_gene_list_count['hidden'] = True
            self._overrides['sample_gene_list_count'] = sample_gene_list_count

        fake_number = "1234567890"
        view_sample_url = reverse('view_sample', kwargs={"sample_id": fake_number}).rstrip(fake_number)
        view_vcf_url = reverse('view_vcf', kwargs={"vcf_id": fake_number}).rstrip(fake_number)
        view_sample_url_prefix = get_url_from_view_path(view_sample_url)
        view_vcf_url_prefix = get_url_from_view_path(view_vcf_url)

        annotation_kwargs = {
            "sample_gene_list_count": Count("samplegenelist", distinct=True),
            "het_hom_count": F("samplestats__het_count") + F("samplestats__hom_count"),
            "sample_url": Func(
                Value(view_sample_url_prefix),
                F("pk"),
                function="CONCAT",
                output_field=CharField(),
            ),
            "vcf_url": Func(
                Value(view_vcf_url_prefix),
                F("vcf_id"),
                function="CONCAT",
                output_field=CharField(),
            ),
        }
        queryset = queryset.annotate(**annotation_kwargs)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, pk):
        """ Do async as it may take a few secs to delete """

        sample = Sample.get_for_user(self.user, pk)
        sample.check_can_write(self.user)
        Sample.objects.filter(pk=sample.pk).update(import_status=ImportStatus.MARKED_FOR_DELETION)
        task = remove_soft_deleted_vcfs_task.si()  # @UndefinedVariable
        task.apply_async()


class CohortSampleListGrid(JqGridUserRowConfig):
    model = Sample
    caption = 'Cohort Samples'
    fields = ["id", "name", "vcf__name", "patient__family_code",
              "patient__first_name", "patient__first_name",
              "patient__sex", "patient__date_of_birth"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewSampleLink'},
                          'vcf__name': {'label': 'VCF'},
                          'patient__family_code': {'label': 'Family Code'},
                          'patient__first_name': {'label': 'First Name'},
                          'patient__last_name': {'label': 'Last Name'},
                          'patient__sex': {'label': 'Sex'},
                          'patient__date_of_birth': {'label': 'D.O.B.'}}

    def __init__(self, user, cohort_id, extra_filters=None):
        super().__init__(user)

        if extra_filters is None:
            extra_filters = {}

        cohort = Cohort.get_for_user(user, cohort_id)
        sample_filters = [Q(vcf__genome_build=cohort.genome_build),
                          Q(import_status=ImportStatus.SUCCESS)]
        SHOW_COHORT = "show_cohort"
        EXCLUDE_COHORT = "exclude_cohort"
        cohort_op = extra_filters.get("cohort_op", EXCLUDE_COHORT)
        cohort_q = Q(cohortsample__cohort=cohort)
        if cohort_op == SHOW_COHORT:
            pass
        elif cohort_op == EXCLUDE_COHORT:
            cohort_q = ~cohort_q
        else:
            raise ValueError(f"Unknown cohort_op: '{cohort_op}'")

        sample_filters.append(cohort_q)
        q = reduce(operator.and_, sample_filters)
        queryset = Sample.filter_for_user(user).filter(q).order_by("-pk")
        self.queryset = queryset.values(*self.get_field_names())


class CohortListGrid(JqGridUserRowConfig):
    model = Cohort
    caption = 'Cohorts'
    fields = ["id", "name", "import_status", "user__username", "modified", "sample_count"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "cohort-icon",
                                      "url_name": "view_cohort",
                                      "url_object_column": "id"}},
        "sample_count": {"label": "Sample Count"},
    }

    def __init__(self, user):
        super().__init__(user)
        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = self.model.filter_for_user(user, success_status_only=False)
        if not user_grid_config.show_group_data:
            queryset = queryset.filter(user=user)
        queryset = queryset.filter(vcf__isnull=True)  # Don't show auto-cohorts

        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': "modified",
                                  'sortorder': "desc"})


class TriosListGrid(JqGridUserRowConfig):
    model = Trio
    caption = 'Trios'
    fields = ["id", "name", "user__username", "modified", "mother__sample__name", "mother_affected",
              "father__sample__name", "father_affected", "proband__sample__name"]
    colmodel_overrides = {
        'id': {'formatter': 'linkFormatter',
               'formatter_kwargs': {"icon_css_class": "trio-icon",
                                    "display_column": "name",
                                    "url_name": "view_trio"}},
        "name": {"hidden": True},
        "mother__sample__name": {"label": "Mother"},
        "father__sample__name": {"label": "Father"},
        "proband__sample__name": {"label": "Proband"}
    }

    def __init__(self, user):
        super().__init__(user)
        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = self.model.filter_for_user(user)
        if not user_grid_config.show_group_data:
            queryset = queryset.filter(user=user)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': "pk",
                                  'sortorder': "desc"})


class GenomicIntervalsListGrid(JqGridUserRowConfig):
    model = GenomicIntervalsCollection
    caption = 'Genomic Intervals'
    fields = ["id", "name", "import_status", "genome_build__name", "user__username"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "bed-icon",
                                      "url_name": "view_genomic_intervals",
                                      "url_object_column": "id"}},
        "genome_build__name": {"label": "Genome Build"},
        'user__username': {'label': 'Uploaded by'}
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = get_objects_for_user(user, 'snpdb.view_genomicintervalscollection', accept_global_perms=False)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())


class CustomColumnsCollectionColumns(DatatableConfig[CustomColumnsCollection]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="user__username", label="User", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    def get_initial_queryset(self) -> QuerySet[CustomColumnsCollection]:
        return CustomColumnsCollection.filter_for_user(self.user)


def server_side_format_clingen_allele(row, field):
    if ca_id := row[field]:
        ca_id = ClinGenAllele.format_clingen_allele(ca_id)
    return ca_id


def server_side_format_exon_and_intron(row, field):
    """ MS Excel will turn '8/11' into a date :( """
    if val := row[field]:
        val = val.replace("/", " of ")
    return val


class AbstractVariantGrid(JqGridUserRowConfig):
    model = Variant

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._count = None
        self.queryset_is_sorted = False

    def _get_standard_overrides(self, af_show_in_percent):
        overrides = {
            # Note:     client side formatters should only be used for adding links etc, never conversion of data, such as
            #           unit to percent, as the CSV downloads (w/o JS formatters) won't match the grid.
            'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink', 'sorttype': 'int'},
            'tags_global': {
                'model_field': False, 'queryset_field': False,
                'name': 'tags_global', 'index': 'tags_global',
                'classes': 'no-word-wrap', 'formatter': 'tagsGlobalFormatter', 'sortable': False
            },
            'clinvar__clinvar_variation_id': {'width': 60, 'formatter': 'clinvarLink'},
            'variantallele__allele__clingen_allele__id': {
                'width': 90,
                "server_side_formatter": server_side_format_clingen_allele,
                'formatter': 'formatClinGenAlleleId'
            },
            'variantannotation__cosmic_id': {'width': 130, 'formatter': 'cosmicLink'},
            'variantannotation__cosmic_legacy_id': {'width': 130, 'formatter': 'cosmicLink'},
            'variantannotation__dbsnp_rs_id': {'width': 130, 'formatter': 'formatDBSNP'},
            'variantannotation__pubmed': {'formatter': 'formatPubMed'},
            'variantannotation__gene__geneannotation__hpo_terms': {'formatter': 'formatOntologyTerms'},
            'variantannotation__gene__geneannotation__mondo_terms': {'formatter': 'formatOntologyTerms'},
            'variantannotation__gene__geneannotation__omim_terms': {'formatter': 'formatOntologyTerms'},
            'variantannotation__transcript_version__gene_version__gene_symbol__symbol': {'formatter': 'geneSymbolLink'},
            'variantannotation__overlapping_symbols': {'formatter': 'geneSymbolNewWindowLink'},
            'variantannotation__transcript_version__gene_version__hgnc__omim_ids': {'width': 60,
                                                                                    'formatter': 'omimLink'},
            'variantannotation__gnomad_filtered': {"formatter": "gnomadFilteredFormatter"},
            'variantannotation__exon': {"server_side_formatter": server_side_format_exon_and_intron},
            'variantannotation__intron': {"server_side_formatter": server_side_format_exon_and_intron},
            'variantannotation__mastermind_mmid3': {'formatter': 'formatMasterMindMMID3'},
            'variantannotation__mavedb_urn': {'formatter': 'formatMavedbUrnLinks'},  # formatMavedbUrnLinks
            # There is more server side formatting (Unit -> Percent) added in _get_fields_and_overrides
        }

        if af_show_in_percent:
            # gnomAD etc are all stored as AF in DB - want to show as percentage on grid
            # But need to be able to turn it off to export VCF as AF
            server_side_format_unit_af = get_allele_frequency_formatter(source_in_percent=False,
                                                                        dest_in_percent=af_show_in_percent)
            af_override = {
                # Unit -> Percent
                'variantannotation__af_1kg': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__af_uk10k': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad2_liftover_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_afr_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_amr_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_asj_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_eas_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_fin_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_nfe_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_oth_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_popmax_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__gnomad_sas_af': {'server_side_formatter': server_side_format_unit_af},
                'variantannotation__topmed_af': {'server_side_formatter': server_side_format_unit_af},
            }
            overrides.update(af_override)
        return overrides

    def _get_base_queryset(self) -> QuerySet:
        raise NotImplementedError()

    def _get_permission_user(self):
        return self.user

    def get_queryset(self, request):
        qs = self._get_base_queryset()
        # Annotate so we can use global_variant_zygosity in grid columns
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = self.filter_items(request, qs)  # JQGrid filtering from request
        if q := self._get_q():
            qs = qs.filter(q)

        field_names = self.get_queryset_field_names()
        a_kwargs = self._get_grid_only_annotation_kwargs()
        qs = qs.annotate(**a_kwargs)
        field_names.extend(a_kwargs)
        return qs.values(*field_names)

    def _get_grid_only_annotation_kwargs(self):
        """ Things not used in counts etc - only to display grid """
        user = self._get_permission_user()
        return get_variantgrid_extra_annotate(user)

    def _get_q(self) -> Optional[Q]:
        return None

    def column_in_queryset_fields(self, field):
        colmodel = self.get_override(field)
        return colmodel.get("queryset_field", True)

    def get_queryset_field_names(self):
        field_names = []
        for f in super().get_field_names():
            if self.column_in_queryset_fields(f):
                field_names.append(f)

        return field_names

    def get_count(self, request):
        if self._count is None:
            queryset = self._get_queryset(request)
            self._count = queryset.count()
        return self._count


class TagColorsCollectionColumns(DatatableConfig[TagColorsCollection]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="user__username", label="User", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    def get_initial_queryset(self) -> QuerySet[TagColorsCollection]:
        return TagColorsCollection.filter_for_user(self.user)


class LiftoverRunColumns(DatatableConfig[LiftoverRun]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            # Don't bother with user as always admin/admin bot
            # RichColumn(key="user__username", label="User", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp',
                       orderable=True, default_sort=SortOrder.DESC),
            # Don't bother with modified as always same time as created
            # RichColumn(key='allele_source', orderable=True),
            RichColumn(key='conversion_tool', renderer=self.render_conversion_tool, orderable=True),
            RichColumn(key='num_alleles', label='Num Alleles', orderable=True, css_class="num"),
            RichColumn(key='source_vcf', orderable=True, css_class="formatted-text"),
            RichColumn(key='source_genome_build', label='Source Build', orderable=True),
            RichColumn(key='genome_build', label='Dest Build', orderable=True),
            RichColumn(key="uploadedliftover__uploaded_file__uploadpipeline__status",
                       label='Status', renderer=self.render_import_status, orderable=True),
            RichColumn(key="uploadedliftover__uploaded_file__uploadpipeline__items_processed",
                       label='Processed', orderable=True, css_class="num"),
        ]

    def render_conversion_tool(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if conversion_tool := row['conversion_tool']:
            act = AlleleConversionTool(conversion_tool)
            label = act.label
        return label

    def render_import_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['uploadedliftover__uploaded_file__uploadpipeline__status']:
            processing_status = ProcessingStatus(status)
            label = processing_status.label
        return label

    def get_initial_queryset(self) -> QuerySet[LiftoverRun]:
        qs = LiftoverRun.objects.all()
        return qs.annotate(num_alleles=Count("alleleliftover"))


class AbstractAlleleLiftoverColumns(DatatableConfig[AlleleLiftover]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request
        self.liftover_run = None

        self.rich_columns = [
            RichColumn(key="allele", orderable=True,
                       renderer=self.render_allele, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="status", label="Status", renderer=self.render_status, orderable=True),
            RichColumn(key="data", label="Data", renderer=self.render_data_json, orderable=True),
            RichColumn(key="error", label="Error", renderer=self.render_error_json, orderable=True),
        ]

    def render_allele(self, row: dict[str, Any]) -> JsonDataType:
        data = {}
        if allele_id := row['allele']:
            allele = get_object_or_404(Allele, id=allele_id)
            data = {
                "text": str(allele),
                "url": allele.get_absolute_url(),
            }
        return data

    def render_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        current = ""
        if status := row['status']:
            processing_status = ProcessingStatus(status)
            label = processing_status.label

        has_build_to_icon = {
            True: "✅",
            False: "❌"
        }

        if allele_id := row['allele']:
            allele: Allele
            if allele := Allele.objects.filter(id=allele_id).first():
                has_37 = bool(allele.variant_for_build_optional(GenomeBuild.grch37()))
                has_38 = bool(allele.variant_for_build_optional(GenomeBuild.grch37()))
                label += f" (Current: {has_build_to_icon[has_37]} GRCh37, {has_build_to_icon[has_38]} GRCh38)"

        return label

    def render_data_json(self, row: dict[str, Any]) -> JsonDataType:
        if js := row["data"]:
            if INFO_LIFTOVER_SWAPPED_REF_ALT in str(js):
                return "Swapped Ref/Alt due to SWAP=1"
        if js is None:
            return "-"
        return jsonify_for_js(js, pretty=True)

    def render_error_json(self, row: dict[str, Any]) -> JsonDataType:
        if js := row["error"]:
            if "message" in js and len(js.keys()) == 1:
                return js.get("message")
        if js is None:
            return "-"
        return jsonify_for_js(js, pretty=True)


class LiftoverRunAlleleLiftoverColumns(AbstractAlleleLiftoverColumns):
    def get_initial_queryset(self) -> QuerySet[AlleleLiftover]:
        liftover_run_id = self.get_query_param("liftover_run_id")
        if liftover_run_id is None:
            raise ValueError("liftover_run_id not provided")
        liftover_run = get_object_or_404(LiftoverRun, pk=liftover_run_id)
        return AlleleLiftover.objects.filter(liftover=liftover_run)


class AlleleLiftoverFailureColumns(AbstractAlleleLiftoverColumns):
    def __init__(self, request):
        super().__init__(request)
        ct_column = RichColumn(key="liftover__conversion_tool",
                               label="Conversion Tool", renderer=self.render_conversion_tool,
                               orderable=True)
        self.rich_columns.insert(1, ct_column)

    def render_conversion_tool(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['liftover__conversion_tool']:
            conversion_tool = AlleleConversionTool(status)
            label = conversion_tool.label
        return label

    def get_initial_queryset(self) -> QuerySet[AlleleLiftover]:
        genome_build_name = self.get_query_param("genome_build_name")
        if genome_build_name is None:
            raise ValueError("genome_build_name not provided")
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        qs = AlleleLiftover.objects.filter(liftover__genome_build=genome_build,
                                           status=ProcessingStatus.ERROR)
        qs = qs.annotate(
            max_id=Max('allele__alleleliftover__id')
        ).filter(
            id=F('max_id')
        )
        return qs


class ManualVariantEntryCollectionColumns(DatatableConfig[ManualVariantEntryCollection]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('manual_variant_entry_collection_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('id', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('user__username', orderable=True),
            RichColumn('import_status', orderable=True, renderer=self._render_import_status),
            RichColumn('genome_build', orderable=True),
        ]

    def _render_import_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['import_status']:
            import_status = ImportStatus(status)
            label = import_status.label
        return label

    def get_initial_queryset(self) -> QuerySet[ManualVariantEntryCollection]:
        qs = ManualVariantEntryCollection.objects.all()
        if not self.user.is_staff:
            qs = qs.filter(user=self.user)
        return qs


class SampleColumns(DatatableConfig[Sample]):
    """ This is currently only used on """
    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id",
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="name", label="Name", orderable=True),
            RichColumn(key="vcf__name", label="VCF", orderable=True),
            RichColumn(key="vcf__date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key=OntologyService.OMIM, orderable=True),
            RichColumn(key=OntologyService.HPO, orderable=True),
            RichColumn(key=OntologyService.MONDO, orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Sample]:
        qs = Sample.filter_for_user(self.user)
        sample_patient_ontology_path = f"patient__{PATIENT_ONTOLOGY_TERM_PATH}"
        ontology_path = f"{sample_patient_ontology_path}__name"
        annotation_kwargs = {}
        for ot in [OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO]:
            q_ot = Q(**{f"{sample_patient_ontology_path}__ontology_service": ot})
            annotation_kwargs[ot.label] = StringAgg(ontology_path, '|',
                                                    filter=q_ot, distinct=True,
                                                    output_field=TextField())
        return qs.annotate(**annotation_kwargs)

    def filter_queryset(self, qs: QuerySet[Sample]) -> QuerySet[Sample]:
        filters = []
        ontology_filters = []
        if ontology_terms := self.get_query_param('ontology_term_id'):
            if q:= get_sample_ontology_q(ontology_terms):
                ontology_filters.append(q)

        if ontology_filters:
            q = reduce(operator.or_, ontology_filters)
            filters.append(q)

        if gene_symbol_str := self.get_query_param("gene_symbol"):
            if q := get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str):
                filters.append(q)

        if filters:
            qs = qs.filter(*filters)
        return qs
