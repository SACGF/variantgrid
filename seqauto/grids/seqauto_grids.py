from functools import cached_property
from typing import Optional

import pandas as pd
from django.db.models import QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse

from genes.models import GeneCoverageCanonicalTranscript
from library.django_utils.datatable_dataframe import DataFrameDatatableConfig
from library.pandas_utils import nan_to_none
from library.utils import JsonDataType, pretty_label
from seqauto.models import EnrichmentKit, GoldCoverageSummary, GoldReference, SequencingSample
from seqauto.seqauto_stats import get_sample_enrichment_kits_df, group_enrichment_kits_df
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


class SequencingSamplesColumns(DatatableConfig[SequencingSample]):
    grid_name = "Sequencing Samples"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="sample_sheet__sequencing_run__name", label="SequencingRun", orderable=True,
                       default_sort=SortOrder.ASC),
            RichColumn(key="sample_name", label="Sample Name", orderable=True),
            RichColumn(key="enrichment_kit__name", label="Enrichment Kit", orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[SequencingSample]:
        return SequencingSample.get_current()


class SequencingSamplesHistoricalConfig(DataFrameDatatableConfig):
    """ Enrichment kit counts over time - one row per time period, one column per kit """
    TIME_FRAMES = {"year": ("Sequencing Samples by year", 'Year'),
                   "month_offset": ("Sequencing Samples by month", 'YYMM'),
                   "all": ("All Sequencing Samples", 'Time Period')}

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.time_frame = self.get_query_param("time_frame")
        if self.time_frame not in self.TIME_FRAMES:
            msg = f"time_frame '{self.time_frame}' invalid! Must be one of: {','.join(self.TIME_FRAMES)}"
            raise ValueError(msg)
        self.csv_name, self.index_label = self.TIME_FRAMES[self.time_frame]

    def get_dataframe(self) -> pd.DataFrame:
        sample_enrichment_kits_df = get_sample_enrichment_kits_df()
        if self.time_frame == "all":
            start = sample_enrichment_kits_df["year_month"].min()
            end = sample_enrichment_kits_df["year_month"].max()
            enrichment_kit_counts = sample_enrichment_kits_df["enrichment_kit__name"].value_counts()
            enrichment_kit_counts.name = f"{start}-{end}"
            df = pd.DataFrame(enrichment_kit_counts).T
        else:
            enrichment_kits_over_time, enrichment_kit_labels = group_enrichment_kits_df(sample_enrichment_kits_df,
                                                                                       self.time_frame)
            enrichment_kit_counts_dict = dict(enrichment_kits_over_time)
            df = pd.DataFrame.from_records(enrichment_kit_counts_dict, index=enrichment_kit_labels)

        return df


class GoldCoverageSummaryColumns(DatatableConfig[GoldCoverageSummary]):
    grid_name = "GoldCoverageSummaries"
    server_csv_download = True
    search_box_enabled = True

    # The search box only makes sense over the symbol/transcript text - the rest are floats
    NUMBER_COLUMNS = ["mean", "standard_error", "min_mean", "depth_20x_5th_percentile",
                      "depth_10x_5th_percentile", "depth_mean_5th_percentile", "depth_mean_95th_percentile"]

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="original_gene_symbol", label="Original Symbol", orderable=True,
                       default_sort=SortOrder.ASC),
            RichColumn(key="original_transcript", label="Original Transcript", orderable=True),
            # ForeignKeys, so the search box has to name the field behind them
            RichColumn(key="gene_symbol", label="Matched Symbol", orderable=True,
                       search=["gene_symbol__symbol"],
                       renderer=self._render_gene_symbol, client_renderer='renderOptionalLink'),
            RichColumn(key="transcript", label="Transcript", orderable=True,
                       search=["transcript__identifier"]),
        ]
        for column in self.NUMBER_COLUMNS:
            renderer = self._render_float if column == "standard_error" else None
            self.rich_columns.append(RichColumn(key=column, label=pretty_label(column), orderable=True,
                                                search=False, css_class="num", renderer=renderer))

    @staticmethod
    def _render_float(cell: CellData) -> Optional[float]:
        return nan_to_none(cell.value)  # JSON can't handle NaN

    def _render_gene_symbol(self, cell: CellData) -> Optional[JsonDataType]:
        if gene_symbol := cell.value:
            return {
                "text": gene_symbol,
                "url": reverse("view_enrichment_kit_gene_coverage",
                               kwargs={"enrichment_kit_id": self.gold_reference.enrichment_kit_id,
                                       "gene_symbol": gene_symbol}),
            }
        return None

    @cached_property
    def gold_reference(self) -> GoldReference:
        return get_object_or_404(GoldReference, pk=self.get_query_param("pk"))

    def get_initial_queryset(self) -> QuerySet[GoldCoverageSummary]:
        return GoldCoverageSummary.objects.filter(gold_reference=self.gold_reference)


class EnrichmentKitGeneCoverageColumns(DatatableConfig[GeneCoverageCanonicalTranscript]):
    grid_name = "Enrichment Kit Gene Coverage"
    server_csv_download = True

    SEQUENCING_SAMPLE_PATH = "gene_coverage_collection__qcgenecoverage__qc__bam_file__sequencing_sample"
    SEQUENCING_RUN_PATH = SEQUENCING_SAMPLE_PATH + "__sample_sheet__sequencing_run"
    GOLD_PATH = SEQUENCING_RUN_PATH + "__gold_standard"
    SAMPLE_NAME_PATH = SEQUENCING_SAMPLE_PATH + "__sample_name"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key=self.SEQUENCING_RUN_PATH, name="sequencing_run", label="Sequencing Run", orderable=True,
                       renderer=self._render_sequencing_run, client_renderer='renderOptionalLink'),
            RichColumn(key=self.GOLD_PATH, name="gold", label="Gold", orderable=True,
                       client_renderer='TableFormat.boolean.bind(null, null)'),
            RichColumn(key=self.SAMPLE_NAME_PATH, name="sample_name", label="Sample", orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key="min", label="Min", orderable=True, css_class="num"),
            RichColumn(key="mean", label="Mean", orderable=True, css_class="num"),
            RichColumn(key="std_dev", label="Std Dev", orderable=True, css_class="num"),
            RichColumn(key="percent_1x", label="% 1x", orderable=True, css_class="num"),
            RichColumn(key="percent_10x", label="% 10x", orderable=True, css_class="num"),
            RichColumn(key="percent_20x", label="% 20x", orderable=True, css_class="num"),
            RichColumn(key="sensitivity", label="Sensitivity", orderable=True, css_class="num"),
        ]

    @staticmethod
    def _render_sequencing_run(cell: CellData) -> Optional[JsonDataType]:
        if sequencing_run := cell.value:
            return {"text": sequencing_run,
                    "url": reverse("view_sequencing_run", kwargs={"sequencing_run_id": sequencing_run})}
        return None

    def get_initial_queryset(self) -> QuerySet[GeneCoverageCanonicalTranscript]:
        enrichment_kit = get_object_or_404(EnrichmentKit, pk=self.get_query_param("enrichment_kit_id"))
        return GeneCoverageCanonicalTranscript.filter_for_kit_and_gene_symbol(
            enrichment_kit, self.get_query_param("genome_build"), self.get_query_param("gene_symbol"))

    def filter_queryset(self, qs: QuerySet[GeneCoverageCanonicalTranscript]) -> QuerySet[GeneCoverageCanonicalTranscript]:
        if self.get_query_param("gold_only") == "true":
            qs = qs.filter(**{self.GOLD_PATH: True})
        return qs
