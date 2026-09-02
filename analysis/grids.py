import operator
from collections import defaultdict
from collections.abc import Callable, Iterator
from functools import reduce
from typing import Any, Optional

import pandas as pd
from auditlog.models import LogEntry
from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.db.models import F, Max, Q, QuerySet, StringAgg, Value
from django.db.models.functions import Substr
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls.base import reverse

from analysis.models import (
    Analysis,
    AnalysisNode,
    AnalysisTemplate,
    AnalysisType,
    Candidate,
    CandidateSearchRun,
    CandidateSearchType,
    CandidateStatus,
    GroupOperation,
    NodeStatus,
)
from analysis.models.models_karyomapping import KaryomappingAnalysis
from analysis.models.nodes.analysis_node import NodeColumnSummaryCacheCollection
from analysis.models.nodes.node_counts import (
    get_extra_filters_count,
    get_node_extra_filters_q,
    is_extra_filter,
)
from analysis.variant_tag_operations import VARIANT_TAG_CLASSIFIED
from analysis.views.analysis_permissions import get_node_subclass_or_404
from annotation.models import HumanProteinAtlasAnnotation
from classification.models import Classification
from genes.grids import GeneListGenesColumns
from genes.models import HGNC, GeneList
from library.django_utils.datatable_dataframe import DataFrameDatatableConfig
from library.unit_percent import get_allele_frequency_formatter
from library.utils import (
    JsonDataType,
    JsonObjType,
    iter_fixed_chunks,
    sha256sum_str,
)
from ontology.grids import AbstractOntologyGenesConfig
from ontology.models import GeneDiseaseClassification, OntologyTermRelation, OntologyVersion
from patients.models_enums import Zygosity
from snpdb.grid_columns.custom_columns import get_variantgrid_extra_annotate
from snpdb.grid_columns.grid_sample_columns import (
    SAMPLE_COMPOSITE_COLUMNS,
    SAMPLE_SORT_KEY_LABELS,
    get_available_format_columns,
    get_variantgrid_zygosity_annotation_kwargs,
)
from snpdb.grids import AbstractVariantGrid
from snpdb.models import (
    CohortGenotype,
    ProcessingStatus,
    Sample,
    UserGridConfig,
    VariantGridColumn,
    VCFFilter,
)
from snpdb.models.models_genome import GenomeBuild
from snpdb.views.datatable_view import CellData, DatatableConfig, NullOrder, RichColumn, SortOrder


class VariantGrid(AbstractVariantGrid):
    grid_name = 'VariantGrid'
    GENOTYPE_COLUMNS_MISSING_VALUE = "."
    SOURCE_COLUMN = "vcf_source"
    # Sorting a packed genotype column orders on an alias annotated for that column alone
    GENOTYPE_SORT_ALIAS_PREFIX = "cohort_genotype_sort_"
    analysis_tags = True
    # Always deferred: whether to fetch rows at all is the page's call (the row count placeholder,
    # and the grid tab being hidden) - node_data_grid.html triggers the first load
    defer_loading = True
    # The fields/operations go down for the FilterNode editor's builder, but the grid itself gets no
    # "Filter grid..." button - FilterNode is how you filter an analysis
    filter_builder_toolbar = False

    def __init__(self, request: HttpRequest, node, extra_filters=None, af_show_in_percent=None):
        self.node = node
        self.genome_build = node.analysis.genome_build
        self.annotation_version = node.analysis.annotation_version
        self.cohorts, self.visibility = node.get_cohorts_and_sample_visibility()
        self.extra_filters = extra_filters
        # column name -> the expression its sort alias is annotated with (@see ordering)
        self._genotype_sort_funcs: dict[str, Any] = {}

        # The rows this view shows - an extra filter (eg clinvar, a tag selection) narrows the node
        self.extra_filters_count = get_extra_filters_count(node, extra_filters)

        super().__init__(request, af_show_in_percent=af_show_in_percent)

    def _get_standard_overrides(self, af_show_in_percent: bool) -> dict[str, dict]:
        overrides = super()._get_standard_overrides(af_show_in_percent)
        # This analysis' own tags - drawn from format_items_iterator/the row, not a queryset field
        overrides['tags'] = {
            'model_field': False, 'queryset_field': False,
            'css_class': 'no-word-wrap', 'orderable': False,
            'client_renderer': 'VariantGridFormat.tags',
        }
        return overrides

    def _get_permission_user(self):
        """ We use the analysis user so it's consistent between users """
        return self.node.analysis.user

    def _get_base_queryset(self) -> QuerySet:
        return self.node.get_queryset()

    def _get_custom_columns_collection(self):
        """ The analysis' columns, so every user sees the same grid """
        return self.node.analysis.custom_columns_collection

    def post_data(self) -> JsonObjType:
        """ The per-request state the page sends as its ajax params on every row request """
        post_data = dict(self.node.get_grid_post_data())
        post_data["node_id"] = self.node.pk
        post_data["version_id"] = self.node.version
        custom_columns_collection = self.node.analysis.custom_columns_collection
        post_data['ccc_id'] = custom_columns_collection.pk
        post_data['ccc_version_id'] = custom_columns_collection.version_id
        post_data["extra_filters"] = self.extra_filters

        sample_ids = self.node.get_sample_ids_with_genotype()
        if sample_ids:
            samples_str = ''.join([str(s) for s in sample_ids])
            post_data['zygosity_samples_hash'] = sha256sum_str(samples_str)
        return post_data

    def _grid_row_count(self) -> Optional[int]:
        """ Current view's row count, or None if unknown """
        if is_extra_filter(self.extra_filters):
            return self.extra_filters_count  # None where we'd have to count the filtered rows
        return self.node.count

    def sorting_disabled(self) -> bool:
        """ Disable sorting on large nodes - sorting by joined/unindexed columns full-sorts the whole
            result set before LIMIT, blowing the statement_timeout (see issue #1651) """
        count = self._grid_row_count()
        return count is None or count >= settings.ANALYSIS_GRID_SORT_MAX_ROWS

    def get_extra(self) -> JsonObjType:
        extra = super().get_extra()
        # Node state the renderers need, once per grid rather than repeated in every column
        extra["analysisNode"] = {"visible": self.node.visible}
        extra["sortingDisabled"] = self.sorting_disabled()
        return extra

    def _get_q(self) -> Optional[Q]:
        return get_node_extra_filters_q(self.node, self.extra_filters)

    def _get_grid_only_annotation_kwargs(self):
        """ Things not used in counts etc - only to display grid """
        annotation_kwargs = get_variantgrid_extra_annotate(self.user, exclude_analysis=self.node.analysis)
        cohorts, _visibility = self.node.get_cohorts_and_sample_visibility()
        common_variants = self.node._has_common_variants()
        try:
            annotation_gnomad_version = self.node.analysis.annotation_version.variant_annotation_version.gnomad
        except AttributeError:
            annotation_gnomad_version = None
        annotation_kwargs.update(get_variantgrid_zygosity_annotation_kwargs(cohorts, common_variants,
                                                                            annotation_gnomad_version=annotation_gnomad_version))
        return annotation_kwargs

    def known_count(self, qs) -> Optional[int]:
        """ The node load pipeline has already counted this node version, so nothing has to count the
            annotated grid queryset again """
        if self.filter_rules_supplied:
            return None  # column filters narrow the rows the stored count was taken over
        count = self._grid_row_count()
        if self.node.count_is_deterministic:
            return count
        # A live-source node's data moves under the stored count. Small nodes are where people curate,
        # and an exact recount at that size is cheap - above it the count stays advisory and the last
        # page may come up short, which beats re-counting millions of rows on every page request
        if count is not None and count <= settings.ANALYSIS_NODE_STORE_ID_SIZE_MAX:
            return None
        return count

    def column(self, name: str) -> RichColumn:
        for rc in self.enabled_columns:
            if rc.name == name:
                return rc
        msg = f"{name} not found in grid columns"
        raise PermissionDenied(msg)

    def ordering(self, qs: QuerySet) -> QuerySet:
        if self.sorting_disabled():
            # Ignore any requested sort - order by -pk only (indexed scan + LIMIT). See issue #1651
            return qs.order_by("-pk")
        for rich_column, _desc in self.requested_ordering():
            # A packed genotype column has no field to sort on - unpack this sample's value into one
            if sort_func := self._genotype_sort_funcs.get(rich_column.name):
                qs = qs.annotate(**{rich_column.sort_keys[0]: sort_func})
        return super().ordering(qs)

    def _get_rich_columns(self) -> list[RichColumn]:
        rich_columns, sample_cols_pos = self._build_variant_grid_columns()

        def insert_columns(columns: list[RichColumn], new_columns: list[RichColumn]) -> list[RichColumn]:
            # Put extra columns after sample (they are all usually to do with sample/vcf etc info)
            new_columns = [rc for rc in new_columns if rc not in columns]
            if sample_cols_pos:
                return columns[:sample_cols_pos] + new_columns + columns[sample_cols_pos:]
            return columns + new_columns

        if extra_columns := self.node.get_extra_columns():
            rich_columns = insert_columns(rich_columns, extra_columns)

        if self.cohorts:
            sample_formatter = None
            if grid_sample_label_template := self.node.analysis.grid_sample_label_template:
                sample_formatter = Sample._get_sample_formatter_func(grid_sample_label_template)

            sample_columns = self._get_grid_genotype_columns(sample_formatter)
            if len(self.cohorts) > 1:
                # Worth a column only where rows can come from more than one VCF
                sample_columns = [self._get_source_column()] + sample_columns
            rich_columns = insert_columns(rich_columns, sample_columns)

        if self.sorting_disabled():
            for rc in rich_columns:
                rc.orderable = False
        elif default_sort_by_column := self.node.analysis.default_sort_by_column:
            # Only set an initial sort column when sorting is allowed (below the row limit) - otherwise
            # the first grid load would request that sort and blow the statement_timeout
            for rc in rich_columns:
                if rc.name == default_sort_by_column.variant_column and rc.orderable:
                    rc.default_sort = SortOrder.ASC
                    break
        return rich_columns

    @staticmethod
    def _get_sample_column_renderer(cohort, sample: Sample, packed_data_replace: dict,
                                    column, i: int, af_show_in_percent: bool):
        """ A function to capture loop variable """
        packed_column = cohort.cohort_genotype_collection.get_packed_column_alias(column)

        def unpack(cell: CellData):
            val = cell[packed_column][i]
            return packed_data_replace.get(val, val)

        renderer = unpack
        if column == "samples_filters" and cohort.vcf:
            filter_formatter = VCFFilter.get_formatter(sample.vcf)  # Per VCF - look up once, not per row

            def sample_filters_renderer(cell: CellData):
                """ Need to unpack then switch filters """
                # Sample Filters can be "."
                val = unpack(cell)
                if val is None:
                    return '.'
                # empty string ('') is PASS
                return filter_formatter(CellData(all_data={packed_column: val}, key=packed_column))

            renderer = sample_filters_renderer
        elif column == "samples_allele_frequency":
            renderer = get_allele_frequency_formatter(source_in_percent=sample.vcf.allele_frequency_percent,
                                                      dest_in_percent=af_show_in_percent,
                                                      get_data_func=unpack,
                                                      missing_value=VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE)
        return renderer

    def _get_sample_cohort_index(self) -> dict:
        """ sample -> (cohort, index into that cohort's packed genotype columns) """
        sample_cohort_index = {}
        for cohort in self.cohorts:
            for cohort_sample in cohort.get_cohort_samples():  # orders by sort_order
                sample = cohort_sample.sample
                if self.visibility.get(sample) and sample not in sample_cohort_index:
                    cohort_index = cohort_sample.cohort_genotype_packed_field_index
                    sample_cohort_index[sample] = (cohort, cohort_index)
        return sample_cohort_index

    def _get_source_column(self) -> RichColumn:
        """ Which VCF a row came from. The pk driven filter of a grouping node carries no provenance,
            but the grid doesn't display from the filter - each row already arrives carrying, per VCF,
            either real packed genotype data or the missing value placeholder. """
        vcf_packed_columns = []
        for sample, (cohort, cohort_index) in self._get_sample_cohort_index().items():
            packed_column = cohort.cohort_genotype_collection.get_packed_column_alias("samples_zygosity")
            vcf_packed_columns.append((str(sample.vcf), packed_column, cohort_index))

        def source_renderer(cell: CellData) -> JsonDataType:
            vcf_names = []
            for vcf_name, packed_column, cohort_index in vcf_packed_columns:
                packed_data = cell.get(packed_column)
                if packed_data and packed_data[cohort_index] != VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE:
                    if vcf_name not in vcf_names:
                        vcf_names.append(vcf_name)
            return ", ".join(vcf_names)

        return RichColumn(key=None, name=VariantGrid.SOURCE_COLUMN, label="Source", width=120,
                          orderable=False,  # Derived from packed data, there's nothing to sort on
                          renderer=source_renderer, csv_rendered=True, include_in_csv=True,
                          extra_columns=sorted({pc for _vcf, pc, _i in vcf_packed_columns}))

    def _get_grid_genotype_columns(self, sample_formatter: Optional[Callable] = None) -> list[RichColumn]:
        available_format_columns = get_available_format_columns(self.cohorts)
        sample_columns = {
            'samples_zygosity': ('Zygosity', '%(sample)s %(label)s', 140),
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
        packed_data_replace.update(dict.fromkeys(MISSING_VALUES, VariantGrid.GENOTYPE_COLUMNS_MISSING_VALUE))

        # We now have separate aliases for packed data, so each cohort handled separately
        rich_columns = []
        for sample, (cohort, cohort_index) in self._get_sample_cohort_index().items():
            for column, (column_label, label_format, width) in sample_columns.items():
                if not available_format_columns[column]:
                    continue
                name = f"sample_{sample.pk}_{column}"
                sample_formatted_str = None
                if sample_formatter:
                    try:
                        sample_formatted_str = sample_formatter(sample)
                    except Exception:
                        pass
                if sample_formatted_str is None or len(sample_formatted_str) == 0:
                    sample_formatted_str = str(sample.name)

                label = label_format % {"sample": sample_formatted_str, "label": column_label}
                renderer = self._get_sample_column_renderer(cohort, sample, packed_data_replace, column,
                                                            cohort_index, self.af_show_in_percent)
                cgc = cohort.cohort_genotype_collection
                self._genotype_sort_funcs[name] = self._genotype_sort_func(cgc, column, sample.pk)
                kwargs = {
                    "key": None,
                    "name": name,
                    "label": label,
                    "width": width,
                    "renderer": renderer,
                    "csv_rendered": True,
                    "include_in_csv": True,
                    "extra_columns": [cgc.get_packed_column_alias(column)],
                    "sort_keys": [self.GENOTYPE_SORT_ALIAS_PREFIX + name],
                    "orderable": True,
                    "null_order": NullOrder.FIRST_ON_ASC,
                }
                if column == 'samples_zygosity':
                    kwargs.update({
                        "client_renderer": 'VariantGridFormat.sampleZygosity',
                        "client_renderer_kwargs": {"samplePrefix": f"sample_{sample.pk}_"},
                        "sort_menu": [
                            {"label": label, "column": f"sample_{sample.pk}_{c}"}
                            for c, label in SAMPLE_SORT_KEY_LABELS.items()
                            if available_format_columns[c]
                        ],
                    })
                elif column in SAMPLE_COMPOSITE_COLUMNS:
                    kwargs["visible"] = False
                rich_columns.append(RichColumn(**kwargs))
        return rich_columns

    @staticmethod
    def _genotype_sort_func(cgc, column: str, sample_id: int):
        """ This sample's value out of the cohort's packed genotype column, as something sortable """
        sql_index = cgc.get_sql_index_for_sample_id(sample_id)
        is_array, _ = CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE[column]
        if is_array:
            # Django index transforms are 0-based
            # https://docs.djangoproject.com/en/5.0/ref/contrib/postgres/fields/#index-transforms
            return F(f"{cgc.cohortgenotype_alias}__{column}__{sql_index - 1}")
        # Is string...
        return Substr(f"{cgc.cohortgenotype_alias}__{column}", sql_index, length=1)


class ExportVariantGrid(VariantGrid):
    """ This is for exporting into VCF/CSV - ie not using any paging """

    EXPORT_PK_BATCH_SIZE = 10000

    def iter_export_rows(self, qs: QuerySet) -> Iterator[dict]:
        """ Take the PKs a contig at a time off the node queryset - that has no annotation joins or
            get_variantgrid_extra_annotate subqueries, so it stays inside the contig/position index -
            then run the full annotated grid queryset against batches of those PKs.

            Every query is bounded by the variant PK index and the annotation joins run once per row at
            any node size. A contig with no variants costs one cheap index probe and no annotated query. """
        value_columns = self.value_columns()
        node_qs = self.node.get_queryset()
        for contig in self.node.analysis.genome_build.standard_contigs:
            contig_pks_qs = node_qs.filter(locus__contig=contig).order_by("locus__position", "pk")
            # A node queryset can fan out over a multi-valued join, so de-dupe (keeping order) to stop
            # a repeated PK straddling a batch boundary and being exported twice
            contig_pks = list(dict.fromkeys(contig_pks_qs.values_list("pk", flat=True)))
            for batch in iter_fixed_chunks(contig_pks, self.EXPORT_PK_BATCH_SIZE):
                batch_qs = qs.filter(pk__in=batch).order_by("locus__position", "pk")
                yield from self.render_export_rows(batch_qs.values(*value_columns).iterator())

    def ordering(self, qs: QuerySet) -> QuerySet:
        """ Export order is set by iter_export_rows (genome build contig, then position), so any
            requested sort is ignored - applying it here would be both wrong and expensive """
        return qs


class AnalysesListColumns(DatatableConfig[Analysis]):
    server_csv_download = True
    # The unfiltered count is over the tag aggregate and the lock join, and only feeds the
    # "(filtered from N total)" text
    count_unfiltered = False

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.scroll_x = True

        self.genome_builds = list(GenomeBuild.builds_with_annotation())
        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, extra_columns=["analysislock__locked"],
                       renderer=self._render_analysis, client_renderer='renderAnalysisLink'),
            RichColumn(key="name", label="Name", orderable=True),
            RichColumn(key="created", label="Created", orderable=True,
                       client_renderer='TableFormat.timestamp'),
            RichColumn(key="modified", label="Modified", orderable=True, default_sort=SortOrder.DESC,
                       client_renderer='TableFormat.timestamp'),
            RichColumn(key="genome_build__name", label="Genome Build", orderable=True,
                       enabled=len(self.genome_builds) > 1),
            RichColumn(key="analysis_type", label="Type", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(AnalysisType.choices)),
            RichColumn(key="description", label="Description", orderable=True),
            RichColumn(key="user__username", label="Created by", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key="tags", label="Tags", client_renderer='renderAnalysisTags'),
            RichColumn(key="id", name="actions", label="", orderable=False,
                       renderer=self._render_actions, client_renderer='renderRowActions'),
        ]

    def _render_analysis(self, cell: CellData) -> JsonDataType:
        analysis_id = cell.value
        return {
            "text": analysis_id,
            "url": reverse("analysis", kwargs={"analysis_id": analysis_id}),
            "locked": bool(cell["analysislock__locked"]),
        }

    def _render_actions(self, cell: CellData) -> JsonDataType:
        analysis_id = cell.value
        actions = [
            {"label": "Settings", "icon": "fas fa-cog", "css_class": "dt-analysis-settings",
             "url": reverse("analysis_settings", kwargs={"analysis_id": analysis_id})},
            {"label": "Clone", "icon": "fas fa-copy", "css_class": "dt-clone-row",
             "url": reverse("clone_analysis", kwargs={"analysis_id": analysis_id})},
        ]
        if delete_url := self.render_delete(cell):
            actions.append({"label": "Delete", "icon": "fas fa-trash", "css_class": "dt-delete-row text-danger",
                            "url": delete_url})
        return actions

    def get_initial_queryset(self) -> QuerySet[Analysis]:
        qs = Analysis.filter_for_user(self.user)
        qs = qs.filter(genome_build__in=self.genome_builds)
        qs = qs.filter(visible=True, template_type__isnull=True)  # Hide templates
        q_last_lock = Q(analysislock=F("last_lock")) | Q(analysislock__isnull=True)
        qs = qs.annotate(last_lock=Max("analysislock__pk")).filter(q_last_lock)
        return qs.annotate(tags=StringAgg("varianttag__tag", delimiter=Value('|')))

    def filter_queryset(self, qs: QuerySet[Analysis]) -> QuerySet[Analysis]:
        user_grid_config = UserGridConfig.get(self.user, 'Analyses')
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs


class AnalysisTemplatesColumns(DatatableConfig[AnalysisTemplate]):
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True, extra_columns=["analysis__id"],
                       renderer=self._render_analysis, client_renderer='renderOptionalLink'),
            RichColumn(key="created", label="Created", orderable=True,
                       client_renderer='TableFormat.timestamp'),
            RichColumn(key="modified", label="Modified", orderable=True, default_sort=SortOrder.DESC,
                       client_renderer='TableFormat.timestamp'),
            RichColumn(key="analysis__genome_build__name", label="Genome Build", orderable=True),
            RichColumn(key="analysis__description", label="Description", orderable=True),
            RichColumn(key="user__username", label="Created by", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key="latest_version", label="Latest version", orderable=True),
            RichColumn(key="active_version", label="Active version", orderable=True),
            RichColumn(key="id", name="clone", label="", orderable=False,
                       renderer=self._render_clone, client_renderer='renderCloneRow'),
            RichColumn(key="id", name="delete", label="", orderable=False,
                       renderer=self.render_delete, client_renderer='TableFormat.deleteRow'),
        ]

    @staticmethod
    def _render_analysis(cell: CellData) -> JsonDataType:
        return {
            "text": cell.value,
            "url": reverse("analysis", kwargs={"analysis_id": cell["analysis__id"]}),
        }

    @staticmethod
    def _render_clone(cell: CellData) -> JsonDataType:
        return reverse("analysis_template_clone", kwargs={"pk": cell.value})

    def get_initial_queryset(self) -> QuerySet[AnalysisTemplate]:
        qs = AnalysisTemplate.filter_for_user(self.user)
        q_active = Q(analysistemplateversion__active=True)
        return qs.annotate(latest_version=Max("analysistemplateversion__version"),
                           active_version=Max("analysistemplateversion__version", filter=q_active))

    def filter_queryset(self, qs: QuerySet[AnalysisTemplate]) -> QuerySet[AnalysisTemplate]:
        user_grid_config = UserGridConfig.get(self.user, 'Analysis Templates')
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs


class NodeColumnSummaryConfig(DataFrameDatatableConfig):
    """ Value counts for one variant grid column across a node's variants """
    index_visible = False
    default_sort_column = "Counts"
    default_sort_order = SortOrder.DESC
    csv_name = "node_column_summary"

    LABELS_COLUMN = "labels"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.node = get_node_subclass_or_404(self.user, self.get_query_param("node_id"),
                                             version=self.get_query_param("node_version"))
        self.extra_filters = self.get_query_param("extra_filters")
        self.variant_column = self.get_query_param("variant_column")
        grid = VariantGrid(request, self.node, self.extra_filters)
        rich_column = grid.column(self.variant_column)
        self.grid_column_name = rich_column.label
        self.renderer = rich_column.renderer

    def get_column_label(self, column_name: Any) -> str:
        if column_name == self.LABELS_COLUMN:
            return self.grid_column_name
        return str(column_name)

    def get_column_client_renderer(self, column_name: Any) -> Optional[str]:
        if column_name == self.LABELS_COLUMN:
            return "renderColumnSummaryFilterChildLink"
        return "TableFormat.number"

    def get_extra(self) -> JsonObjType:
        """ The labels renderer only links where the column can become a FilterNode """
        return {
            "variantColumn": self.variant_column,
            "createFilterChildLinks": VariantGridColumn.objects.filter(variant_column=self.variant_column).exists(),
        }

    def get_dataframe(self) -> pd.DataFrame:
        counts = NodeColumnSummaryCacheCollection.get_counts_for_node(self.node, self.variant_column,
                                                                     self.extra_filters)
        if self.renderer:
            # The same server side rendering the grid cell got (choices expanded, AF as percent)
            labels = {field: self.renderer(CellData(all_data={self.variant_column: field},
                                                    key=self.variant_column))
                      for field in counts}
        else:
            labels = {field: field for field in counts}

        counts_series = pd.Series(counts)
        df = pd.DataFrame({self.LABELS_COLUMN: labels, "Counts": counts_series})
        total = counts_series.sum()
        if total != 0:
            df['Percent'] = 100.0 * df['Counts'] / total
        else:
            df['Percent'] = 0.0
        return df.sort_values("Counts", ascending=False)


class NodeOntologyGenesConfig(AbstractOntologyGenesConfig):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.node = get_node_subclass_or_404(self.user, self.get_query_param("node_id"),
                                             version=self.get_query_param("version"))

    def _get_ontology_version(self) -> OntologyVersion:
        return self.node.analysis.annotation_version.ontology_version

    def _get_ontology_term_ids(self):
        return self.node.get_ontology_term_ids()


class NodeGeneDiseaseClassificationGenesConfig(DataFrameDatatableConfig):
    """ Gene symbols in a MOI node, with each source term's mode of inheritance summary """
    index_label = "Gene Symbol"
    csv_name = "gene_disease_classification_genes"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.node = get_node_subclass_or_404(self.user, self.get_query_param("node_id"),
                                             version=self.get_query_param("version"))

    def _get_ontology_term_relations(self) -> list[OntologyTermRelation]:
        return self.node.get_gene_disease_relations()

    def get_dataframe(self) -> pd.DataFrame:
        gene_data = defaultdict(dict)
        valid_classifications = list(reversed(GeneDiseaseClassification.labels))
        for otr in self._get_ontology_term_relations():
            moi_classifications = otr.get_gene_disease_moi_classifications()
            gene_symbol = otr.dest_term.name
            summary = ", ".join(otr.get_moi_summary(moi_classifications, valid_classifications))
            gene_data[gene_symbol][str(otr.source_term)] = summary

        df = pd.DataFrame.from_dict(gene_data, orient='index')
        return df.sort_index()


class AnalysisNodeIssuesColumns(DatatableConfig[AnalysisNode]):
    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='analysis__id', visible=False),
            RichColumn(key='analysis__name', label='Analysis', orderable=True,
                       renderer=self.render_node_link,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='status', orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(NodeStatus.choices)),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key='errors', orderable=False),
        ]

    def render_node_link(self, row: CellData) -> JsonDataType:
        analysis_id = row['analysis__id']
        node_id = row['id']
        url = reverse('analysis_node', kwargs={'analysis_id': analysis_id, 'active_node_id': node_id})
        return {"text": row.value, "url": url}

    def get_initial_queryset(self) -> QuerySet[AnalysisNode]:
        return AnalysisNode.objects.filter(status=NodeStatus.ERROR)


class KaryomappingAnalysesColumns(DatatableConfig[KaryomappingAnalysis]):
    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key='trio__cohort__genome_build__name', label='Genome Build', orderable=True),
            RichColumn(key="user__username", label='Created by', orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key='trio__name', label='Trio', orderable=True),
            RichColumn(key='trio__proband__sample__name', label='Proband', orderable=True),
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[KaryomappingAnalysis]:
        return KaryomappingAnalysis.filter_for_user(self.user)

    def filter_queryset(self, qs: QuerySet[KaryomappingAnalysis]) -> QuerySet[KaryomappingAnalysis]:
        user_grid_config = UserGridConfig.get(self.user, 'Karomapping Analyses')
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs


class NodeTissueExpressionGenesColumns(DatatableConfig[HumanProteinAtlasAnnotation]):
    download_csv_button_enabled = True

    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='gene_symbol', label='Gene Symbol', orderable=True),
            RichColumn(key='gene', label='Gene', orderable=True),
            RichColumn(key='value', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[HumanProteinAtlasAnnotation]:
        node_id = self.get_query_param('node_id')
        version = self.get_query_param('version')
        node = get_node_subclass_or_404(self.user, node_id, version=version)
        return node.get_hpa_qs()


class NodeTissueUniProtGenesColumns(DatatableConfig[HGNC]):
    download_csv_button_enabled = True

    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='gene_symbol__symbol', label='Gene Symbol', orderable=True),
            RichColumn(key='uniprot__accession', label='UniProt Accession', orderable=True),
            RichColumn(key='uniprot__tissue_specificity', label='Tissue Specificity', orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[HGNC]:
        node_id = self.get_query_param('node_id')
        version = self.get_query_param('version')
        node = get_node_subclass_or_404(self.user, node_id, version=version)
        filters = []
        for word in node.text_tissue.split():
            filters.append(Q(uniprot__tissue_specificity__icontains=word))
        q = GroupOperation.reduce(filters, node.group_operation)
        return HGNC.objects.filter(q)


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

    if content_type_model == "varianttag":
        if additional_data.get("operation") == VARIANT_TAG_CLASSIFIED:
            classification_id = additional_data["classification_id"]
            url = Classification.get_url_for_pk(classification_id)
            return f"Retired - classified as <a href='{url}'>{classification_id}</a>"

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
            RichColumn(key="actor__username", orderable=True,
                       extra_columns=["actor__id"], renderer=self.render_user),
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
            RichColumn(key="user__username", label="User", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
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
                data["url"] = reverse("classify_candidate", args=[row["id"]])
                data["text"] = "📑"
                data["title"] = "Classify sample"

        return data

    @staticmethod
    def render_clinvar(row: dict[str, Any]) -> JsonDataType:
        data = {}
        if evidence := row.get("evidence"):
            data = evidence.get("ClinVar diff", {})
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
            RichColumn(key="user__username", label='User', orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
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
