import operator
import time
from collections import defaultdict
from functools import reduce
from typing import Optional, Any, Callable

import pandas as pd
from auditlog.models import LogEntry
from django.conf import settings
from django.contrib.postgres.aggregates import StringAgg
from django.core.exceptions import PermissionDenied
from django.db.models import Max, F, Q, QuerySet
from django.db.models.functions import Substr
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django.utils.functional import SimpleLazyObject

from analysis.models import Analysis, AnalysisNode, NodeCount, NodeStatus, AnalysisTemplate, GroupOperation, \
    CandidateSearchRun, CandidateSearchType, Candidate, CandidateStatus, AnalysisType
from analysis.models.models_karyomapping import KaryomappingAnalysis
from analysis.models.nodes.analysis_node import get_extra_filters_q, NodeColumnSummaryCacheCollection
from analysis.views.analysis_permissions import get_node_subclass_or_404
from annotation.models import HumanProteinAtlasAnnotation
from genes.grids import GeneListGenesColumns
from genes.models import HGNC, GeneList
from library.jqgrid.jqgrid_sql import get_overrides
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.pandas_jqgrid import DataFrameJqGrid
from library.unit_percent import get_allele_frequency_formatter
from library.utils import update_dict_of_dict_values, JsonDataType, sha256sum_str
from ontology.grids import AbstractOntologyGenesGrid
from ontology.models import OntologyTermRelation, GeneDiseaseClassification, OntologyVersion
from patients.models_enums import Zygosity
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position, \
    get_variantgrid_extra_annotate
from snpdb.grid_columns.grid_sample_columns import get_available_format_columns, \
    get_variantgrid_zygosity_annotation_kwargs
from snpdb.grids import AbstractVariantGrid
from snpdb.models import VariantGridColumn, UserGridConfig, VCFFilter, Sample, CohortGenotype, ProcessingStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class VariantGrid(AbstractVariantGrid):
    caption = 'VariantGrid'
    GENOTYPE_COLUMNS_MISSING_VALUE = "."
    colmodel_overrides = {
        'tags': {'classes': 'no-word-wrap', 'formatter': 'tagsFormatter', 'sortable': False},
    }

    def __init__(self, user, node, extra_filters=None, af_show_in_percent=None):
        if af_show_in_percent is None:
            af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT

        self.cohorts, self.visibility = node.get_cohorts_and_sample_visibility()
        self.fields, override = self._get_fields_and_overrides(node, af_show_in_percent)
        super().__init__(user)  # Need to call init after setting fields

        self.url = SimpleLazyObject(lambda: reverse("node_grid_handler", kwargs={"analysis_id": node.analysis_id}))
        self.extra_config.update(node.get_extra_grid_config())
        default_sort_by_column = node.analysis.default_sort_by_column
        if default_sort_by_column:
            self.extra_config['sortname'] = default_sort_by_column.variant_column

        update_dict_of_dict_values(self._overrides, override)

        self.node = node
        self.name = node.name

        try:
            node_count = NodeCount.load_for_node(self.node, extra_filters)
        except NodeCount.DoesNotExist:
            node_count = None
        self.node_count = node_count
        self._set_post_data(node, extra_filters)

    def _get_permission_user(self):
        """ We use the analysis user so it's consistent between users """
        return self.node.analysis.user

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
            post_data['zygosity_samples_hash'] = sha256sum_str(samples_str)
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

    def _get_grid_only_annotation_kwargs(self):
        """ Things not used in counts etc - only to display grid """
        annotation_kwargs = get_variantgrid_extra_annotate(self.user, exclude_analysis=self.node.analysis)
        cohorts, _visibility = self.node.get_cohorts_and_sample_visibility()
        common_variants = self.node._has_common_variants()
        annotation_kwargs.update(get_variantgrid_zygosity_annotation_kwargs(cohorts, common_variants))
        return annotation_kwargs

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

    def _get_fields_and_overrides(self, node: AnalysisNode, af_show_in_percent: bool) -> tuple[list, dict]:
        ccc = node.analysis.custom_columns_collection
        annotation_version = node.analysis.annotation_version
        fields, overrides, sample_cols_pos = get_custom_column_fields_override_and_sample_position(ccc,
                                                                                                   annotation_version,
                                                                                                   analysis_tags=True)
        # Put extra columns after sample (they are all usually to do with sample/vcf etc info)
        if extra_columns := node.get_extra_columns():
            if sample_cols_pos:
                fields = fields[:sample_cols_pos] + extra_columns + fields[sample_cols_pos:]
            else:
                fields.extend(extra_columns)

        update_dict_of_dict_values(overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(overrides, node.get_extra_colmodel_overrides())
        if self.cohorts:
            sample_formatter = None
            if grid_sample_label_template := node.analysis.grid_sample_label_template:
                sample_formatter = Sample._get_sample_formatter_func(grid_sample_label_template)

            sample_columns, sample_overrides = VariantGrid.get_grid_genotype_columns_and_overrides(self.cohorts, self.visibility,
                                                                                                   af_show_in_percent, sample_formatter)
            if sample_cols_pos:
                fields = fields[:sample_cols_pos] + sample_columns + fields[sample_cols_pos:]
            else:
                fields.extend(sample_columns)

            update_dict_of_dict_values(overrides, sample_overrides)
        return fields, overrides

    @staticmethod
    def _get_sample_columns_server_side_formatter(cohort, sample: Sample, packed_data_replace: dict,
                                                  column, i: int, af_show_in_percent: bool):
        """ A function to capture loop variable """

        def packed_data_formatter(row, _field):
            cgc = cohort.cohort_genotype_collection
            packed_column = cgc.get_packed_column_alias(column)
            packed_data = row[packed_column]
            val = packed_data[i]
            return packed_data_replace.get(val, val)

        server_side_formatter = packed_data_formatter
        if column == "samples_filters" and cohort.vcf:
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
    def get_grid_genotype_columns_and_overrides(cohorts, visibility,
                                                af_show_in_percent: bool, sample_formatter: Optional[Callable] = None):
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
        MISSING_VALUES = [CohortGenotype.MISSING_NUMBER_VALUE, CohortGenotype.MISSING_FT_VALUE, -2147483648]
        packed_data_replace.update({mv: VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE for mv in MISSING_VALUES})

        # We now have separate aliases for packed data, so each cohort handled separately
        sample_cohort_index = {}
        for cohort in cohorts:
            for cohort_sample in cohort.get_cohort_samples():  # orders by sort_order
                sample = cohort_sample.sample
                if visibility.get(sample) and sample not in sample_cohort_index:
                    cohort_index = cohort_sample.cohort_genotype_packed_field_index
                    sample_cohort_index[sample] = (cohort, cohort_index)

        column_names = []
        column_data = []
        for sample, (cohort, cohort_index) in sample_cohort_index.items():
            for column, (column_label, label_format, width) in sample_columns.items():
                if not available_format_columns[column]:
                    continue
                column_names.append(f"sample_{sample.pk}_{column}")
                sample_formatted_str = None
                if sample_formatter:
                    try:
                        sample_formatted_str = sample_formatter(sample)
                    except:
                        pass
                if sample_formatted_str is None or len(sample_formatted_str) == 0:
                    sample_formatted_str = str(sample.name)

                label = label_format % {"sample": sample_formatted_str, "label": column_label}
                server_side_formatter = VariantGrid._get_sample_columns_server_side_formatter(cohort, sample,
                                                                                              packed_data_replace,
                                                                                              column, cohort_index,
                                                                                              af_show_in_percent)
                cgc = cohort.cohort_genotype_collection
                sql_index = cgc.get_sql_index_for_sample_id(sample.pk)

                col_data_dict = {
                    "label": label,
                    "width": width,
                    "server_side_formatter": server_side_formatter,
                    # Index is what is passed back to server side for sorting - we'll pack the info here
                    "index": ":".join([cgc.cohortgenotype_alias, str(sql_index), column]),
                }
                column_data.append(col_data_dict)

        overrides = get_overrides(column_names, column_data, model_field=False, queryset_field=False)
        return column_names, overrides

    def _sort_items(self, items, sidx, sord):
        """ Special case to handle sort by CohortGenotype packed fields """
        if sidx is not None:
            # For special fields, we pack sorting info into the 'index' which doesn't map to a field
            # looks like 'cohortgenotype_134:1:samples_zygosity'
            if ":" in sidx:
                sort_alias = "cohort_genotype_sample_sort_alias"
                cohortgenotype_alias, sql_index, column = sidx.split(":")
                sql_index = int(sql_index)
                is_array, _ = CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE[column]
                if is_array:
                    # Django index transforms are 0-based
                    # https://docs.djangoproject.com/en/5.0/ref/contrib/postgres/fields/#index-transforms
                    django_index = sql_index-1
                    sort_func = F(f"{cohortgenotype_alias}__{column}__{django_index}")
                else:
                    # Is string...
                    sort_func = Substr(f"{cohortgenotype_alias}__{column}", sql_index, length=1)

                items = items.annotate(**{sort_alias: sort_func})
                sidx = sort_alias

        return super()._sort_items(items, sidx, sord)


class ExportVariantGrid(VariantGrid):
    """ This is for exporting into VCF/CSV - ie not using any paging """

    def __init__(self, *args, **kwargs):
        self.order_by = kwargs.pop("order_by", None)
        super().__init__(*args, **kwargs)

    def sort_items(self, request, items):
        if self.order_by:
            items = items.order_by(*self.order_by)
        return items

    @staticmethod
    def _iter_time(items):
        start = time.time()
        yield from items
        end = time.time()
        print(f"Download took {end-start} seconds")

    @staticmethod
    def _iter_by_contigs(genome_build, items):
        for contig in genome_build.standard_contigs:
            # print(f"getting {contig}")
            contig_items = items.filter(locus__contig=contig).iterator()
            yield from contig_items

    def paginate_items(self, request, items):
        # This is the step after queryset is sorted
        # We want everything, but will retrieve contig at a time to reduce DB query
        genome_build = self.node.analysis.genome_build
        items = self._iter_by_contigs(genome_build, items)
        # items = self._iter_time(items)
        return None, None, items


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
        "tags": {"label": "Tags", "model_field": False,
                 "name": "tags", "formatter": "tagsFormatter"},  # This formatter counts multiple tags
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

    def _get_ontology_term_relations(self) -> list[OntologyTermRelation]:
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


class NodeGeneListGenesColumns(GeneListGenesColumns):
    """ This shows internals of gene lists, using permission of analysis instead of request.user for gene list
        It's assumed that if someone added the gene list via auto-complete then they give people who can see
        that analysis permissions (according to analysis not gene list) """

    def _get_gene_annotation_releases(self) -> list['GeneAnnotationRelease']:
        analysis_id = self.get_query_param("analysis_id")
        analysis = Analysis.get_for_user(self.user, pk=analysis_id)
        return [analysis.gene_annotation_release]

    def _get_gene_list(self):
        gene_list_id = self.get_query_param("gene_list_id")
        node_id = self.get_query_param("node_id")
        version = self.get_query_param("version")

        gene_list = get_object_or_404(GeneList, pk=gene_list_id)
        node = get_node_subclass_or_404(self.user, node_id, version=version)
        gene_list_in_node = False
        for gl in node.get_gene_lists():
            if gl == gene_list:
                gene_list_in_node = True
                break
        if not gene_list_in_node:
            raise PermissionDenied(f"GeneList {gene_list_id} not part of GeneListNode gene lists")
        return gene_list


def get_analysis_log_entry_summary(action, content_type_model, changes, additional_data) -> Optional[str]:
    """ If it can come up with a summary, return something, otherwise nothing, and the
        caller can do it themselves """
    if content_type_model == "analysisedge":
        parent = additional_data["parent"]
        child = additional_data["child"]
        parent_desc = f"(id={parent['id']},type:{parent['content_type']['model']})"
        if parent_rep := parent['object_repr']:
            parent_desc += f" '{parent_rep}'"
        child_desc = f"(id={child['id']},type:{child['content_type']['model']})"
        if child_rep := child['object_repr']:
            child_desc += f" '{child_rep}'"

        if action == LogEntry.Action.CREATE:
            op = "Attached"
            desc = "to"
        elif action == LogEntry.Action.DELETE:
            op = "Detached"
            desc = "from"
        else:
            raise ValueError("Don't know how to handle analysisedge UPDATE")
        return f"{op} Parent: {parent_desc} {desc} Child: {child_desc}"

    if action == LogEntry.Action.CREATE:
        return "Created"
    elif action == LogEntry.Action.DELETE:
        return "Deleted"
    elif action == LogEntry.Action.UPDATE:
        # Just a move operation
        changes_set = set(changes.keys())
        move_set = {"x", "y"}
        # If only x/y it's a move
        if changes_set | move_set and not (changes_set - move_set):
            x = changes.get("x")
            y = changes.get("y")
            if x is not None:
                if y is not None:
                    return f"Moved to {changes['x'][1]},{changes['y'][1]}"
                else:
                    return f"Moved to x={x[1]}"
            else:
                return f"Moved to y={y[1]}"

        is_save = changes_set == {"count", "status", "version", "appearance_version"}
        is_save &= changes.get("status", [None, None])[1] == NodeStatus.DIRTY
        if is_save:
            return "Saved"
    return None


class AnalysisLogEntryColumns(DatatableConfig[LogEntry]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user
        self.analysis = None

        self.expand_client_renderer = "renderExpandAnalysisAuditLogEntry"
        self.rich_columns = [
            RichColumn(key="timestamp", orderable=True, client_renderer='TableFormat.timestamp'),
            RichColumn(key="actor__username", orderable=True),
            RichColumn(key="content_type__model", orderable=True),
            RichColumn(key=None, name="node_id", label="Node ID", renderer=self.render_node_id),
            RichColumn(key="object_repr", label="Object"),
            RichColumn(key="action", label="Action", renderer=self.render_action),
            RichColumn(key=None, name='summary', label='Summary',
                       renderer=self.render_summary, client_renderer='renderAnalysisAuditLogSummary'),
            RichColumn(key="changes", visible=False),
            RichColumn(key="additional_data", visible=False),
        ]

    def render_action(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        action = row['action']
        if action is not None:
            lookup = dict(LogEntry.Action.choices)
            label = lookup[action]
        return label

    def render_node_id(self, row: dict[str, Any]) -> JsonDataType:
        return row["additional_data"].get("node_id")

    def render_summary(self, row: dict[str, Any]) -> JsonDataType:
        action = row['action']
        content_type_model = row["content_type__model"]
        changes = row["changes"]
        additional_data = row["additional_data"]

        summary_json = {}
        if summary_text := get_analysis_log_entry_summary(action, content_type_model, changes, additional_data):
            summary_json["summary_text"] = summary_text
        else:
            summary_json['changes'] = changes
            summary_json['additional_data'] = additional_data
        return summary_json

    def get_initial_queryset(self) -> QuerySet[LogEntry]:
        if node_id := self.get_query_param("node_id"):
            node = get_node_subclass_or_404(self.user, node_id)
            qs = node.log_entry_qs()
        else:
            analysis_id = self.get_query_param("analysis_id")
            if analysis_id is None:
                raise ValueError("'analysis_id' not provided")
            analysis = Analysis.get_for_user(self.user, pk=analysis_id)
            qs = analysis.log_entry_qs()
        return qs


class CandidateSearchRunColumns(DatatableConfig[LogEntry]):
    def __init__(self, request):
        super().__init__(request)
        self.user = request.user
        self.rich_columns = [
            RichColumn('id',
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl', default_sort=SortOrder.DESC),
            RichColumn(key="search_version__search_type", orderable=True, renderer=self.render_search_type),
            RichColumn(key="search_version__code_version", orderable=True),
            RichColumn(key="created", label="Created", orderable=True, client_renderer='TableFormat.timestamp'),
            RichColumn(key="user__username", label="User", orderable=True),
            RichColumn(key="status", label="Status", orderable=True, renderer=self.render_status),
        ]

    def get_initial_queryset(self) -> QuerySet[CandidateSearchRun]:
        qs = CandidateSearchRun.filter_for_user(self.user)
        if search_types := self.get_query_param("search_types"):
            qs = qs.filter(search_version__search_type__in=search_types)
        return qs

    @staticmethod
    def render_search_type(row: dict[str, Any]):
        return CandidateSearchType(row["search_version__search_type"]).label

    @staticmethod
    def render_status(row: dict[str, Any]):
        return ProcessingStatus(row["status"]).label


class CandidateColumns(DatatableConfig[LogEntry]):
    def __init__(self, request):
        super().__init__(request)
        self.user = request.user
        csr_id = self.get_query_param("candidate_search_run_id")
        # Retrieve for Permission check
        self.csr = CandidateSearchRun.get_for_user(self.user, pk=csr_id)

        self.rich_columns = [
            RichColumn('id', visible=False),
            RichColumn(key="status", orderable=True, renderer=self.render_status),
            RichColumn(name="action", label="Action",
                       renderer=self.render_action, client_renderer='candidate_action_renderer'),
            RichColumn(key="variant", label="Variant", orderable=True,
                       renderer=self.render_variant_link, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="notes", orderable=True),
            RichColumn(key="evidence", label="Evidence", visible=False, orderable=True, client_renderer='TableFormat.json'),
            # RichColumn(key="reviewer__username", label="Reviewer", orderable=True),
            # RichColumn(key="reviewer_comment", label="Reviewer Comment", orderable=True),
            RichColumn(
                key='sample_id',
                name='sample_id',
                visible=False,  # Only used to build links
            ),
        ]

        # Show/hide various columns based on search type (as we only use some)
        optional_columns = [
            RichColumn(key="classification__clinical_significance",
                       label="Clin Sig",
                       orderable=True,
                       client_renderer="classification_clinical_significance_renderer"),
            RichColumn(key="classification",
                       label="Classification",
                       orderable=True,
                       extra_columns=[
                           'classification__evidence__c_hgvs__value',
                           'classification__evidence__g_hgvs__value',
                           'classification__condition_resolution__display_text',
                       ],
                       renderer=self.render_classification_summary,
                       client_renderer="classification_summary_renderer"),
            RichColumn(key="sample__name", label="Sample", orderable=True, client_renderer='VCTable.sample'),
            RichColumn(key="analysis", label="Analysis", orderable=True,
                       renderer=self.render_analysis_link, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="annotation_version", label="Annotation Version", orderable=True),
            RichColumn(key="zygosity", label="Zygosity", orderable=True, renderer=self.render_zygosity),
            RichColumn(name="clinvar", label="ClinVar",
                       renderer=self.render_clinvar, client_renderer='clinvar_renderer'),
        ]
        columns = CandidateSearchRun.CANDIDATE_GRID_COLUMNS[self.csr.search_version.search_type]
        for rc in optional_columns:
            if rc.name in columns:
                self.rich_columns.append(rc)

    def get_initial_queryset(self) -> QuerySet[Candidate]:
        qs = Candidate.objects.filter(search_run=self.csr)

        if candidate_status := self.get_query_param("candidate_status"):
            candidates = candidate_status.split(",")
            qs = qs.filter(status__in=candidates)

        if evidence := self.get_query_param("evidence"):
            evidence_keys = evidence.split(",")
            q_list = [Q(**{f"evidence__{ek}__isnull": False}) for ek in evidence_keys]
            if q_list:
                q = reduce(operator.or_, q_list)
                qs = qs.filter(q)

        # Aggregate filters
        for col in ["sample", "analysis", "classification"]:
            if value := self.get_query_param(col):
                qs = qs.filter(**{col: value})
        return qs

    @staticmethod
    def render_status(row: dict[str, Any]):
        return CandidateStatus(row["status"]).label

    @staticmethod
    def render_variant_link(row: dict[str, Any]):
        variant = row["variant"]
        text = variant
        if g_hgvs := row.get("classification__evidence__g_hgvs__value"):
            text = g_hgvs

        return {
            "text": text,
            "url": reverse('view_variant', kwargs={'variant_id': variant}),
        }

    @staticmethod
    def render_analysis_link(row: dict[str, Any]):
        data = {}
        if analysis := row.get("analysis"):
            data = {
                "text": analysis,
                "url": reverse('analysis', kwargs={'analysis_id': analysis}),
            }
        return data

    @staticmethod
    def render_zygosity(row: dict[str, Any]):
        return Zygosity.display(row["zygosity"])

    @staticmethod
    def render_classification_summary(row: dict[str, Any]) -> JsonDataType:
        return {
            "classification": row["classification"],
            'c_hgvs': row.get('classification__evidence__c_hgvs__value'),
            'g_hgvs': row.get('classification__evidence__g_hgvs__value'),
            'classification__condition_resolution__display_text': row.get('classification__condition_resolution__display_text'),
        }

    @staticmethod
    def render_action(row: dict[str, Any]) -> JsonDataType:
        data = {}
        if row["status"] in (CandidateStatus.OPEN, CandidateStatus.HIGHLIGHTED):
            # Only do this if there is a sample
            if row.get("sample_id"):
                data["url"] = reverse("classify_candidate", args=[row["id"]]),
                data["text"] = "📑"
                data["title"] = "Classify sample"

        return data

    @staticmethod
    def render_clinvar(row: dict[str, Any]) -> JsonDataType:
        data = {}
        if evidence := row.get("evidence"):
            data = evidence.get("clinvar", {})
        return data


class AnalysesColumns(DatatableConfig[Analysis]):
    """ For listing Analyses """
    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn('id',
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="name", orderable=True),
            RichColumn(key="user__username", label='User', orderable=True),
            RichColumn(key="analysis_type", label="Type", orderable=True, renderer=self.render_analysis_type),
            RichColumn(key="created", label='Created', orderable=True, client_renderer='TableFormat.timestamp'),
            RichColumn(key="modified", label='Modified', orderable=True, client_renderer='TableFormat.timestamp'),

        ]

    def get_initial_queryset(self) -> QuerySet[Analysis]:
        qs = Analysis.filter_for_user(self.user)
        qs = qs.filter(visible=True, template_type__isnull=True)  # Hide templates

        params = ['analysis_type', 'my_user', 'date_min', 'date_max']
        data = {k: self.get_query_param(k) for k in params}

        if filters := self.get_q_list(self.user, data):
            qs = qs.filter(*filters)

        return qs

    @staticmethod
    def get_q_list(user, params: dict) -> list[Q]:
        """ This has been split off so that reanalysis tasks can use the same filters """
        q_list = []
        if analysis_type := params.get('analysis_type'):
            q_list.append(Q(analysis_type=analysis_type))

        if params.get('my_user'):
            q_list.append(Q(user=user))

        if date_min := params.get('date_min'):
            q_list.append(Q(created__gte=date_min))

        if date_max := params.get('date_max'):
            q_list.append(Q(modified__lte=date_max))
        return q_list

    @staticmethod
    def render_analysis_type(row: dict[str, Any]):
        if analysis_type := row["analysis_type"]:
            analysis_type = AnalysisType(row["analysis_type"]).label
        return analysis_type

