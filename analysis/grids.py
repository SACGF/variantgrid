from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import pandas as pd
from django.conf import settings
from django.contrib.postgres.aggregates import StringAgg
from django.core.exceptions import PermissionDenied
from django.db.models import Max, F, Q, QuerySet
from django.urls.base import reverse
from django.utils.functional import SimpleLazyObject

from analysis.models import Analysis, AnalysisNode, NodeCount, NodeStatus, AnalysisTemplate, GroupOperation
from analysis.models.models_karyomapping import KaryomappingAnalysis
from analysis.models.nodes.analysis_node import get_extra_filters_q, NodeColumnSummaryCacheCollection
from analysis.views.analysis_permissions import get_node_subclass_or_404
from annotation.models import HumanProteinAtlasAnnotation
from genes.models import HGNC
from library.jqgrid.jqgrid_sql import get_overrides
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.pandas_jqgrid import DataFrameJqGrid
from library.unit_percent import get_allele_frequency_formatter
from library.utils import md5sum_str
from ontology.grids import AbstractOntologyGenesGrid
from ontology.models import OntologyTermRelation, GeneDiseaseClassification, OntologyVersion
from patients.models_enums import Zygosity
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position, \
    get_variantgrid_extra_alias_and_select_columns
from snpdb.grid_columns.grid_sample_columns import get_columns_and_sql_parts_for_cohorts, get_available_format_columns
from snpdb.grids import AbstractVariantGrid
from snpdb.models import VariantGridColumn, UserGridConfig, VCFFilter, Sample, CohortGenotype
from snpdb.models.models_genome import GenomeBuild


class VariantGrid(AbstractVariantGrid):
    caption = 'VariantGrid'
    GENOTYPE_COLUMNS_MISSING_VALUE = "."
    colmodel_overrides = {
        'tags': {'classes': 'no-word-wrap', 'formatter': 'tagsFormatter', 'sortable': False},
    }

    def __init__(self, user, node, extra_filters=None, sort_by_contig_and_position=False, af_show_in_percent=None):
        if af_show_in_percent is None:
            af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT

        self.fields, override = self._get_fields_and_overrides(node, af_show_in_percent)
        super().__init__(user)  # Need to call init after setting fields

        self.url = SimpleLazyObject(lambda: reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id}))
        self.sort_by_contig_and_position = sort_by_contig_and_position
        self.extra_config.update(node.get_extra_grid_config())
        default_sort_by_column = node.analysis.default_sort_by_column
        if default_sort_by_column:
            self.extra_config['sortname'] = default_sort_by_column.variant_column

        self.update_overrides(override)

        self.node = node
        self.name = node.name

        try:
            node_count = NodeCount.load_for_node(self.node, extra_filters)
        except NodeCount.DoesNotExist:
            node_count = None
        self.node_count = node_count
        self._set_post_data(node, extra_filters)

    def _get_base_queryset(self) -> QuerySet:
        return self.node.get_queryset()

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

    def _get_q(self) -> Optional[Q]:
        q = None
        if self.node_count:
            analysis = self.node.analysis
            q = get_extra_filters_q(analysis.user, analysis.genome_build, self.node_count.label)
        return q

    def _get_variantgrid_extra_alias_and_select_columns(self):
        return get_variantgrid_extra_alias_and_select_columns(self.user, exclude_analysis=self.node.analysis)

    def _get_new_columns_select_from_where_parts(self, values_queryset) -> Tuple[List[str], str, str, str]:
        cohorts, _visibility = self.node.get_cohorts_and_sample_visibility()

        if cohorts:
            ret = get_columns_and_sql_parts_for_cohorts(values_queryset, cohorts)
            new_columns, select_part, from_part, where_part = ret
        else:
            return super()._get_new_columns_select_from_where_parts(values_queryset)

        return new_columns, select_part, from_part, where_part

    def get_count(self, request):  # pylint: disable=unused-argument
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

    def _get_fields_and_overrides(self, node: AnalysisNode, af_show_in_percent: bool) -> Tuple[List, Dict]:
        ccc = node.analysis.custom_columns_collection
        annotation_version = node.analysis.annotation_version
        fields, overrides, sample_columns_position = get_custom_column_fields_override_and_sample_position(ccc, annotation_version)
        fields.extend(node.get_extra_columns())
        overrides.update(node.get_extra_colmodel_overrides())
        overrides.update(self._get_standard_overrides(af_show_in_percent))
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
              "user__username", "tags", "analysislock__locked"]
    colmodel_overrides = {
        'id': {'formatter': 'analysisLink',
               'formatter_kwargs': {"icon_css_class": "analysis-icon",
                                    "url_name": "analysis"}},
        "name": {"width": 500},
        "genome_build__name": {"label": "Genome Build"},
        "analysis_type": {"label": "Type"},
        "user__username": {'label': 'Created by'},
        "tags": {"label": "Tags", "model_field": False, "formatter": "tagsFormatter"},  # This formatter counts multiple tags
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
        qs = Analysis.filter_for_user(user)
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=user)
        qs = qs.filter(genome_build__in=self.genome_builds)
        qs = qs.filter(visible=True, template_type__isnull=True)  # Hide templates
        q_last_lock = Q(analysislock=F("last_lock")) | Q(analysislock__isnull=True)
        qs = qs.annotate(last_lock=Max("analysislock__pk")).filter(q_last_lock)
        qs = qs.annotate(tags=StringAgg("varianttag__tag", delimiter='|'))
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

    def __init__(self, user, analysis_id, node_id, node_version, extra_filters, variant_column, significant_figures):
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

        if VariantGridColumn.objects.filter(variant_column=variant_column).exists():
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
    fields = ["id", "analysis__id", "analysis__name", "status", "modified", "errors"]
    colmodel_overrides = {
        "id": {"hidden": True},
        'analysis__id': {"hidden": True},
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
    fields = ['id', 'name', 'modified', "trio__cohort__genome_build__name", 'user__username',
              'trio__name', 'trio__proband__sample__name']

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
        queryset = KaryomappingAnalysis.filter_for_user(user)
        if not user_grid_config.show_group_data:
            queryset = queryset.filter(user=user)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'modified', 'sortorder': 'desc'})


class NodeOntologyGenesGrid(AbstractOntologyGenesGrid):
    colmodel_overrides = {
        "ID": {"width": 200},
        "hpo": {"width": 400},
        "omim": {"width": 400},
    }

    def __init__(self, user, analysis_id, node_id, version):
        self.node = get_node_subclass_or_404(user, node_id, version=version)
        super().__init__()

    def _get_ontology_version(self) -> OntologyVersion:
        return self.node.analysis.annotation_version.ontology_version

    def _get_ontology_term_ids(self):
        return self.node.get_ontology_term_ids()


class NodeGeneDiseaseClassificationGenesGrid(DataFrameJqGrid):
    def __init__(self, user, analysis_id, node_id, version):
        super().__init__()
        self.node = get_node_subclass_or_404(user, node_id, version=version)

    def _get_ontology_term_relations(self) -> List[OntologyTermRelation]:
        return self.node.get_gene_disease_relations()

    def get_dataframe(self):
        gene_data = defaultdict(dict)
        valid_classifications = list(reversed(GeneDiseaseClassification.labels))
        columns = {}
        for otr in self._get_ontology_term_relations():
            moi_classifications = otr.get_gene_disease_moi_classifications()
            gene_symbol = otr.dest_term.name
            summary = ", ".join(otr.get_moi_summary(moi_classifications, valid_classifications))
            column = str(otr.source_term)
            columns[column] = True
            gene_data[gene_symbol][column] = summary

        self._overrides.update({column: {"width": 500} for column in columns})

        df = pd.DataFrame.from_dict(gene_data, orient='index')
        return df.sort_index()


class NodeTissueExpressionGenesGrid(JqGridUserRowConfig):
    model = HumanProteinAtlasAnnotation
    caption = 'Tissue Node: Human Protein Atlas'
    fields = ['gene_symbol', 'gene', 'value']

    def __init__(self, user, analysis_id, node_id, version):
        super().__init__(user)
        node = get_node_subclass_or_404(user, node_id, version=version)
        queryset = node.get_hpa_qs()
        self.queryset = queryset.values(*self.get_field_names())


class NodeTissueUniProtTissueSpecificityGenesGrid(JqGridUserRowConfig):
    model = HGNC
    caption = 'Tissue Node: UniProt Tissue Specificity'
    fields = ['gene_symbol__symbol', 'uniprot__accession', 'uniprot__tissue_specificity']

    def __init__(self, user, analysis_id, node_id, version):
        super().__init__(user)
        node = get_node_subclass_or_404(user, node_id, version=version)
        filters = []
        for word in node.text_tissue.split():
            f = Q(uniprot__tissue_specificity__icontains=word)
            filters.append(f)
        q = GroupOperation.reduce(filters, node.group_operation)
        queryset = self.model.objects.filter(q)
        self.queryset = queryset.values(*self.get_field_names())
