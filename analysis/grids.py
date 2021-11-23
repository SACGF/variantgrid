import pandas as pd
from typing import Dict, List, Tuple

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.db.models import Max, F, Q
from django.urls.base import reverse
from django.utils.functional import SimpleLazyObject

from analysis.models import Analysis, AnalysisNode, NodeCount, NodeStatus, AnalysisTemplate
from analysis.models.models_karyomapping import KaryomappingAnalysis
from analysis.models.nodes.analysis_node import get_extra_filters_q, NodeColumnSummaryCacheCollection
from analysis.views.analysis_permissions import get_node_subclass_or_404
from library.database_utils import get_queryset_column_names, get_queryset_select_from_where_parts
from library.jqgrid_sql import JqGridSQL, get_overrides
from library.jqgrid_user_row_config import JqGridUserRowConfig
from library.pandas_jqgrid import DataFrameJqGrid
from library.unit_percent import get_allele_frequency_formatter
from library.utils import md5sum_str
from ontology.grids import AbstractOntologyGenesGrid
from patients.models_enums import Zygosity
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position, \
    get_variantgrid_extra_alias_and_select_columns
from snpdb.grid_columns.grid_sample_columns import get_columns_and_sql_parts_for_cohorts, get_available_format_columns
from snpdb.models import VariantGridColumn, UserGridConfig, VCFFilter, Sample, ClinGenAllele, CohortGenotype
from snpdb.models.models_genome import GenomeBuild


def server_side_format_clingen_allele(row, field):
    if ca_id := row[field]:
        ca_id = ClinGenAllele.format_clingen_allele(ca_id)
    return ca_id


def server_side_format_exon_and_intron(row, field):
    """ MS Excel will turn '8/11' into a date :( """
    if val := row[field]:
        val = val.replace("/", " of ")
    return val


class VariantGrid(JqGridSQL):
    model = AnalysisNode.model
    caption = 'VariantGrid'
    url = SimpleLazyObject(lambda: reverse("node_grid_handler"))
    GENOTYPE_COLUMNS_MISSING_VALUE = "."
    colmodel_overrides = {
        # Note:     client side formatters should only be used for adding links etc, never conversion of data, such as
        #           unit to percent, as the CSV downloads (w/o JS formatters) won't match the grid.
        'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink', 'sorttype': 'int'},
        'tags': {'classes': 'no-word-wrap', 'formatter': 'tagsFormatter', 'sortable': False},
        'tags_global': {'classes': 'no-word-wrap', 'formatter': 'tagsGlobalFormatter', 'sortable': False},
        'clinvar__clinvar_variation_id': {'width': 60, 'formatter': 'clinvarLink'},
        'variantallele__allele__clingen_allele': {'width': 90,
                                                  "server_side_formatter": server_side_format_clingen_allele,
                                                  'formatter': 'formatClinGenAlleleId'},
        'variantannotation__cosmic_id': {'width': 90, 'formatter': 'cosmicLink'},
        'variantannotation__cosmic_legacy_id': {'width': 90, 'formatter': 'cosmicLink'},
        'variantannotation__transcript_version__gene_version__gene_symbol__symbol': {'formatter': 'geneSymbolLink'},
        'variantannotation__overlapping_symbols': {'formatter': 'geneSymbolNewWindowLink'},
        'variantannotation__transcript_version__gene_version__hgnc__omim_ids': {'width': 60, 'formatter': 'omimLink'},
        'variantannotation__gnomad_filtered': {"formatter": "gnomadFilteredFormatter"},
        'variantannotation__exon': {"server_side_formatter": server_side_format_exon_and_intron},
        'variantannotation__intron': {"server_side_formatter": server_side_format_exon_and_intron},
        # There is more server side formatting (Unit -> Percent) added in _get_fields_and_overrides
    }

    def __init__(self, user, node, extra_filters=None, sort_by_contig_and_position=False, af_show_in_percent=None):
        if af_show_in_percent is None:
            af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT

        self.fields, override = self._get_fields_and_overrides(node, af_show_in_percent)
        super().__init__(user)  # Need to call init after setting fields

        self.sort_by_contig_and_position = sort_by_contig_and_position
        self.extra_config.update(node.get_extra_grid_config())
        default_sort_by_column = node.analysis.default_sort_by_column
        if default_sort_by_column:
            self.extra_config['sortname'] = default_sort_by_column.variant_column

        self.update_overrides(override)

        self.node = node
        self.name = node.name
        self.queryset_is_sorted = False

        try:
            node_count = NodeCount.load_for_node(self.node, extra_filters)
        except NodeCount.DoesNotExist:
            node_count = None
        self.node_count = node_count
        self._set_post_data(node, extra_filters)

    def _set_post_data(self, node, extra_filters):
        post_data = self.extra_config.get('postData', {})
        post_data["node_id"] = node.pk
        post_data["version_id"] = node.version
        custom_columns_collection = node.analysis.custom_columns_collection
        post_data['ccc_id'] = custom_columns_collection.pk
        post_data['ccc_version_id'] = custom_columns_collection.version_id
        post_data["extra_filters"] = extra_filters

        sample_ids = node.get_sample_ids()
        if sample_ids:
            samples_str = ''.join([str(s) for s in sample_ids])
            zygosity_samples_hash = md5sum_str(samples_str)
            post_data['zygosity_samples_hash'] = zygosity_samples_hash
        self.extra_config['postData'] = post_data

    def get_colmodels(self, remove_server_side_only=False):
        """ Put 'analysisNode' into every colmodel """
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        global_colmodel = {"analysisNode": {"visible": self.node.visible}}
        for cm in colmodels:
            cm.update(global_colmodel)
        return colmodels

    def _get_queryset(self):
        qs = self.node.get_queryset()
        if self.sort_by_contig_and_position:
            # These 2 are custom_columns.ID_FORMATTER_REQUIRED_FIELDS so won't add extra cols
            qs = qs.order_by("locus__contig__name", "locus__position")
            self.queryset_is_sorted = True
        return qs

    def get_values_queryset(self):
        self.queryset = self._get_queryset()
        if self.node_count:
            analysis = self.node.analysis
            extra_filters_q = get_extra_filters_q(analysis.user, analysis.genome_build, self.node_count.label)
            self.queryset = self.queryset.filter(extra_filters_q)
        self.queryset = self.queryset.values(*self.get_queryset_field_names())
        return self.queryset

    def column_in_queryset_fields(self, field):
        colmodel = self.get_override(field)
        return colmodel.get("queryset_field", True)

    def get_queryset_field_names(self):
        field_names = []
        for f in super().get_field_names():
            if self.column_in_queryset_fields(f):
                field_names.append(f)

        return field_names

    def get_sql_params_and_columns(self, request):
        values_queryset = self.get_values_queryset()
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

        sql, column_names = self._get_zygosity_sql_and_columns(values_queryset, order_by)
        return sql, [], column_names, is_sorted

    def _get_zygosity_sql_and_columns(self, values_queryset, order_by: str):
        cohorts, _visibility = self.node.get_cohorts_and_sample_visibility()

        if cohorts:
            ret = get_columns_and_sql_parts_for_cohorts(values_queryset, cohorts)
            new_columns, select_part, from_part, where_part = ret
        else:
            new_columns = []
            select_part, from_part, where_part = get_queryset_select_from_where_parts(values_queryset)

        extra_column_selects = []
        for alias, select in get_variantgrid_extra_alias_and_select_columns(self.user,
                                                                            exclude_analysis=self.node.analysis):
            extra_column_selects.append(f'({select}) as "{alias}"')
            new_columns.append(alias)

        if order_by:  # SELECT DISTINCT, ORDER BY expressions must appear in select list
            extra_column_selects.append(order_by)
        select_part += ",\n" + ",\n".join(extra_column_selects)

        sql = '\n'.join([select_part, from_part, where_part])
        column_names = get_queryset_column_names(values_queryset, new_columns)
        return sql, column_names

    def get_count(self):
        """ Used by paginator, set from stored value so that we don't
            have to make another SQL query """
        if self.node_count:
            count = self.node_count.count
        elif self.node.count is not None:
            count = self.node.count
        else:
            class_name = self.node.get_class_name()
            node_details_list = [f"pk: {self.node.pk}",
                                 f"version: {self.node.version}",
                                 f"status: {self.node.status}"]
            node_details = ", ".join(node_details_list)
            msg = f"Node {class_name} ({node_details}) does not have count set"
            raise ValueError(msg)
        return count

    def get_column_colmodel(self, column_name):
        for cm in self.get_colmodels():
            if column_name == cm['name']:
                return cm

        msg = f"{column_name} not found in grid column model"
        raise PermissionDenied(msg)

    @staticmethod
    def _get_fields_and_overrides(node: AnalysisNode, af_show_in_percent: bool) -> Tuple[List, Dict]:
        ccc = node.analysis.custom_columns_collection
        fields, overrides, sample_columns_position = get_custom_column_fields_override_and_sample_position(ccc)
        fields.extend(node.get_extra_columns())
        overrides.update(node.get_extra_colmodel_overrides())
        if af_show_in_percent:
            # gnomAD etc are all stored as AF in DB
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

        cohorts, visibility = node.get_cohorts_and_sample_visibility()
        if cohorts:
            sample_columns, sample_overrides = VariantGrid.get_grid_genotype_columns_and_overrides(cohorts, visibility, af_show_in_percent)
            if sample_columns_position:
                fields = fields[:sample_columns_position] + sample_columns + fields[sample_columns_position:]
            else:
                fields.extend(sample_columns)
            overrides.update(sample_overrides)
        return fields, overrides

    @staticmethod
    def _get_sample_columns_server_side_formatter(sample: Sample, packed_data_replace: Dict,
                                                  column, i: int, af_show_in_percent: bool):
        """ A function to capture loop variable """

        def packed_data_formatter(row, _field):
            packed_data = row[f"packed_{column}"]
            val = packed_data[i]
            return packed_data_replace.get(val, val)

        server_side_formatter = packed_data_formatter
        if column == "samples_filters":
            def sample_filters_formatter(row, field):
                """ Need to unpack then switch filters """
                # Sample Filters can be "."
                val = packed_data_formatter(row, field)
                if val is None:
                    return '.'
                # empty string ('') is PASS
                filter_formatter = VCFFilter.get_formatter(sample.vcf)
                row[field] = val
                return filter_formatter(row, field)

            server_side_formatter = sample_filters_formatter
        elif column == "samples_allele_frequency":
            server_side_formatter = get_allele_frequency_formatter(source_in_percent=sample.vcf.allele_frequency_percent,
                                                                   dest_in_percent=af_show_in_percent,
                                                                   get_data_func=packed_data_formatter,
                                                                   missing_value=VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE)
        return server_side_formatter

    @staticmethod
    def get_grid_genotype_columns_and_overrides(cohorts, visibility, af_show_in_percent: bool):
        available_format_columns = get_available_format_columns(cohorts)
        sample_columns = {
            'samples_zygosity': ('Zygosity', '%(sample)s %(label)s', 55),
            'samples_allele_depth': ('AD', '%(label)s %(sample)s', 25),
            'samples_allele_frequency': ('AF', '%(label)s %(sample)s', 30),
            'samples_read_depth': ('DP', '%(label)s %(sample)s', 25),
            'samples_genotype_quality': ('GQ', '%(label)s %(sample)s', 25),
            'samples_phred_likelihood': ('PL', '%(label)s %(sample)s', 25),
            'samples_filters': ('FT', '%(label)s %(sample)s', 100),
        }
        packed_data_replace = dict(Zygosity.CHOICES)
        # Some legacy data (Missing data in FreeBayes before PythonKnownVariantsImporter v12) has -2147483647 for
        # empty values (what CyVCF2 returns using format()) @see https://github.com/SACGF/variantgrid/issues/59
        MISSING_VALUES = [CohortGenotype.MISSING_NUMBER_VALUE, -2147483648]
        packed_data_replace.update({mv: VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE for mv in MISSING_VALUES})

        # Record the 1st cohort a sample appears in, and the concatenated cohorts index
        sample_cohort_cat_cohorts_index = {}
        cohorts_offset = 0
        for cohort in cohorts:
            for cohort_sample in cohort.get_cohort_samples():  # orders by sort_order
                sample = cohort_sample.sample
                if visibility.get(sample) and sample not in sample_cohort_cat_cohorts_index:
                    cc_index = cohorts_offset + cohort_sample.cohort_genotype_packed_field_index
                    sample_cohort_cat_cohorts_index[sample] = (cohort, cc_index)
            cohorts_offset += cohort.cohort_genotype_collection.num_samples

        column_names = []
        column_data = []
        for sample, (cohort, cc_index) in sample_cohort_cat_cohorts_index.items():
            for column, (column_label, label_format, width) in sample_columns.items():
                if not available_format_columns[column]:
                    continue
                column_names.append(f"{sample.pk}_{column}")
                label = label_format % {"sample": sample.name, "label": column_label}
                server_side_formatter = VariantGrid._get_sample_columns_server_side_formatter(sample,
                                                                                              packed_data_replace,
                                                                                              column, cc_index,
                                                                                              af_show_in_percent)
                col_data_dict = {
                    "label": label,
                    "width": width,
                    "server_side_formatter": server_side_formatter,
                    "order_by": cohort.get_sample_column_order_by(sample, column)
                }
                column_data.append(col_data_dict)

        overrides = get_overrides(column_names, column_data, model_field=False, queryset_field=False)
        return column_names, overrides


class AnalysesGrid(JqGridUserRowConfig):
    model = Analysis
    caption = 'Analyses'
    fields = ["id", "name", "created", "modified", "genome_build__name", "analysis_type", "description",
              "user__username", "analysislock__locked"]
    colmodel_overrides = {
        'id': {'formatter': 'analysisLink',
               'formatter_kwargs': {"icon_css_class": "analysis-icon",
                                    "url_name": "analysis"}},
        "name": {"width": 500},
        "genome_build__name": {"label": "Genome Build"},
        "analysis_type": {"label": "Type"},
        "user__username": {'label': 'Created by'},
        "analysislock__locked": {"hidden": True},
    }

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)
        fields = self.get_field_names()

        self.genome_builds = list(GenomeBuild.builds_with_annotation())
        if len(self.genome_builds) == 1:  # No need to show
            genome_build_colmodel = self._overrides.get('genome_build', {})
            genome_build_colmodel['hidden'] = True
            self._overrides['genome_build'] = genome_build_colmodel
        user_grid_config = UserGridConfig.get(user, self.caption)
        if user_grid_config.show_group_data:
            qs = Analysis.filter_for_user(user)
        else:
            qs = Analysis.objects.filter(user=user)
        qs = qs.filter(genome_build__in=self.genome_builds)
        qs = qs.filter(visible=True, template_type__isnull=True)  # Hide templates
        q_last_lock = Q(analysislock=F("last_lock")) | Q(analysislock__isnull=True)
        qs = qs.annotate(last_lock=Max("analysislock__pk")).filter(q_last_lock)
        self.queryset = qs.values(*fields)
        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})


class AnalysisTemplatesGrid(JqGridUserRowConfig):
    model = AnalysisTemplate
    caption = 'Analysis Templates'
    fields = ["id", "analysis__id", "name", "created", "modified",
              "analysis__genome_build__name", "analysis__description", "user__username"]

    colmodel_overrides = {
        "id": {"hidden": True},  # Need an ID row so we can delete
        'analysis__id': {'formatter': 'analysisLink',
                         'formatter_kwargs': {"icon_css_class": "analysis-icon",
                                              "url_name": "analysis"}},
        "name": {"width": 500},
        "analysis": {"hidden": True},
        "modified": {'label': 'Modified'},
        "user__username": {'label': 'Created by'},
    }

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)

        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = AnalysisTemplate.filter_for_user(user)
        if not user_grid_config.show_group_data:
            queryset = queryset.filter(user=user)
        queryset = queryset.annotate(latest_version=Max("analysistemplateversion__version")).values("latest_version")
        fields = self.get_field_names() + ["latest_version"]
        self.queryset = queryset.values(*fields)
        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        extra = {'index': 'latest_version', 'name': 'latest_version', 'label': 'Latest version', 'sorttype': 'int'}
        colmodels.append(extra)
        return colmodels

    def delete_row(self, pk):
        analysis_template = AnalysisTemplate.get_for_user(self.user, pk)
        analysis_template.check_can_write(self.user)
        analysis_template.delete_or_soft_delete()


class NodeColumnSummaryGrid(DataFrameJqGrid):
    colmodel_overrides = {
        "ID": {"hidden": True},
        "labels": {'formatter': 'createFilterChildLink'},
        "Percent": {"formatter": "number"},
    }

    def __init__(self, user, node_id, node_version, extra_filters, variant_column, significant_figures):
        super().__init__()

        self.node = get_node_subclass_or_404(user, node_id, version=node_version)
        grid = VariantGrid(user, self.node, extra_filters)
        cm = grid.get_column_colmodel(variant_column)
        grid_column_name = cm["label"]
        field_formatters = grid.get_field_formatters()
        self.formatter = field_formatters.get(variant_column)
        self.extra_filters = extra_filters
        self.variant_column = variant_column
        self.significant_figures = significant_figures
        self.grid_column_name = grid_column_name
        self._overrides["labels"]["label"] = grid_column_name
        self.extra_config.update({'sortname': 'Counts',
                                  'sortorder': 'desc'})

        # Only want to
        variantgrid_column = VariantGridColumn.objects.filter(variant_column=variant_column).first()
        if variantgrid_column:
            self.extra_config["create_filter_child_links"] = True

    def get_dataframe(self):
        counts = NodeColumnSummaryCacheCollection.get_counts_for_node(self.node, self.variant_column, self.extra_filters)
        if self.formatter:
            labels = {}
            for field in counts:
                fake_row = {"field": field}
                labels[field] = self.formatter(fake_row, "field")
        else:
            labels = {field: field for field in counts}

        counts_series = pd.Series(counts)
        df = pd.DataFrame({"labels": labels, "Counts": counts_series})
        total = counts_series.sum()
        if total != 0:
            df['Percent'] = 100.0 * df['Counts'] / total
        else:
            df['Percent'] = 0.0
        return df.sort_values("Counts", ascending=False)


class AnalysisNodeIssuesGrid(JqGridUserRowConfig):
    model = AnalysisNode
    caption = 'Analysis Node Issues'
    fields = ["id", "analysis_id", "analysis__name", "status", "modified", "errors"]
    colmodel_overrides = {
        "id": {"hidden": True},
        'analysis_id': {"hidden": True},
        'analysis__name': {"width": 400,
                           'formatter': 'analysisNodeLink',
                           'formatter_kwargs': {"icon_css_class": "analysis-icon",
                                                "url_name": "analysis_node",
                                                "url_object_column": "hack_kwargs"}},
        "modified": {'label': 'Modified'},
    }

    def __init__(self, **kwargs):
        user = kwargs.get("user")
        super().__init__(user)

        queryset = self.model.objects.filter(status=NodeStatus.ERROR)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})


class KaromappingAnalysesGrid(JqGridUserRowConfig):
    model = KaryomappingAnalysis
    caption = 'Karomapping Analyses'
    fields = ['id', 'name', 'modified', "trio__cohort__genome_build", 'user__username', 'trio__name', 'trio__proband__sample__name']

    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {"width": 400,
                 'formatter': 'linkFormatter',
                 'formatter_kwargs': {"url_name": "view_karyomapping_analysis",
                                      "url_object_column": "id"}},
        "trio__cohort__genome_build": {"label": "Genome Build"},
        "user__username": {'label': 'Created by'},
        "trio__name": {"label": "Trio"},
        "trio__proband__sample__name": {"label": "Proband"}
    }

    def __init__(self, user):
        super().__init__(user)

        user_grid_config = UserGridConfig.get(user, self.caption)
        if user_grid_config.show_group_data:
            queryset = KaryomappingAnalysis.filter_for_user(user)
        else:
            queryset = KaryomappingAnalysis.objects.filter(user=user)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'modified', 'sortorder': 'desc'})


class NodeOntologyGenesGrid(AbstractOntologyGenesGrid):
    colmodel_overrides = {
        "ID": {"width": 200},
        "hpo": {"width": 400},
        "omim": {"width": 400},
    }

    def __init__(self, user, node_id, version):
        self.node = get_node_subclass_or_404(user, node_id, version=version)
        super().__init__()

    def _get_ontology_term_ids(self):
        return self.node.get_ontology_term_ids()
