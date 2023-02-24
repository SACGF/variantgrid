import operator
from functools import reduce
from typing import Optional, List, Tuple

from django.conf import settings
from django.db.models import F, QuerySet
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from guardian.shortcuts import get_objects_for_user

from library.jqgrid.jqgrid_sql import JqGridSQL
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.unit_percent import get_allele_frequency_formatter
from library.utils import calculate_age
from library.utils.database_utils import get_queryset_column_names, \
    get_queryset_select_from_where_parts, queryset_to_sql
from snpdb.grid_columns.custom_columns import get_variantgrid_extra_annotate
from snpdb.models import VCF, Cohort, Sample, ImportStatus, \
    GenomicIntervalsCollection, CustomColumnsCollection, Variant, Trio, UserGridConfig, GenomeBuild, ClinGenAllele, \
    VariantZygosityCountCollection, TagColorsCollection
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs, remove_soft_deleted_vcfs_task
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class VCFListGrid(JqGridUserRowConfig):
    model = VCF
    caption = 'VCFs'
    fields = ["id", "name", "date", "import_status", "genome_build__name", "user__username", "source",
              "uploadedvcf__uploaded_file__import_source", "genotype_samples", "project__name", "cohort__import_status",
              "uploadedvcf__vcf_importer__name", 'uploadedvcf__vcf_importer__version']
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'width': 550,
                 'formatter': 'viewVCFLink',
                 'formatter_kwargs': {"icon_css_class": "vcf-icon",
                                      "url_name": "view_vcf",
                                      "url_object_column": "id"}},
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

        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, vcf_id):
        """ Do async as it may be slow """
        soft_delete_vcfs(self.user, vcf_id)


# TODO: Merge this an cohort grid below into 1
class SamplesListGrid(JqGridUserRowConfig):
    model = Sample
    caption = 'Samples'
    fields = ["id", "name", "het_hom_count", "vcf__date", "import_status", "vcf__genome_build__name", "variants_type",
              "vcf__user__username", "vcf__source", "vcf__name", "vcf__project__name", "vcf__uploadedvcf__uploaded_file__import_source",
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
        "specimen__reference_id": {'label': 'Specimen'}
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

        annotation_kwargs = {
            "sample_gene_list_count": Count("samplegenelist", distinct=True),
            "het_hom_count": F("samplestats__het_count") + F("samplestats__hom_count"),
        }
        queryset = queryset.annotate(**annotation_kwargs)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, sample_id):
        """ Do async as it may take a few secs to delete """

        sample = Sample.get_for_user(self.user, sample_id)
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

    def get_initial_queryset(self) -> QuerySet[TagColorsCollection]:
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


class AbstractVariantGrid(JqGridSQL):
    model = Variant

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._count = None
        self.queryset_is_sorted = False
        self.sort_by_contig_and_position = False

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

    def _get_queryset(self, request):
        qs = self._get_base_queryset()
        # Annotate so we can use global_variant_zygosity in grid columns
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = self.filter_items(request, qs)  # JQGrid filtering from request
        if q := self._get_q():
            qs = qs.filter(q)

        if self.sort_by_contig_and_position:
            # These 2 are custom_columns.ID_FORMATTER_REQUIRED_FIELDS so won't add extra cols
            qs = qs.order_by("locus__contig__name", "locus__position")
            self.queryset_is_sorted = True
        return qs

    def get_values_queryset(self, request, field_names: List = None):
        queryset = self._get_queryset(request)
        if field_names is None:
            field_names = self.get_queryset_field_names()
            a_kwargs = self._get_grid_only_annotation_kwargs()
            queryset = queryset.annotate(**a_kwargs)
            field_names.extend(a_kwargs)
        return queryset.values(*field_names)

    def get_sql_params_and_columns(self, request):
        values_queryset = self.get_values_queryset(request)
        is_sorted = self.queryset_is_sorted
        order_by = None
        if not is_sorted:
            # If sort column is normal, go through normal sort_items path.
            # If special, modify SQL below
            sidx = request.GET.get('sidx', 'id')
            if self.column_in_queryset_fields(sidx):
                values_queryset = self.sort_items(request, values_queryset)
                is_sorted = True

            override = self.get_override(sidx)
            order_by = override.get("order_by")

        sql, column_names = self._get_sql_and_columns(values_queryset, order_by)
        return sql, [], column_names, is_sorted

    def _get_grid_only_annotation_kwargs(self):
        """ Things not used in counts etc - only to display grid """
        return get_variantgrid_extra_annotate(self.user)

    def _get_sql_and_columns(self, values_queryset, order_by: str):
        new_columns, select_part, from_part, where_part = self._get_new_columns_select_from_where_parts(values_queryset)

        extra_column_selects = []
        if order_by:  # SELECT DISTINCT, ORDER BY expressions must appear in select list
            extra_column_selects.append(order_by)
        if extra_column_selects:
            select_part += ",\n" + ",\n".join(extra_column_selects)

        sql = '\n'.join([select_part, from_part, where_part])
        column_names = get_queryset_column_names(values_queryset, new_columns)
        return sql, column_names

    def _get_new_columns_select_from_where_parts(self, values_queryset) -> Tuple[List[str], str, str, str]:
        select_part, from_part, where_part = get_queryset_select_from_where_parts(values_queryset)
        return [], select_part, from_part, where_part

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
