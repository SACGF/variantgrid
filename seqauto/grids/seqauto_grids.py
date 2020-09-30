from django.shortcuts import get_object_or_404
from django.urls.base import reverse

from genes.models import GeneCoverageCanonicalTranscript
from library.django_utils.jqgrid_view import JQGridViewOp
from library.jqgrid_user_row_config import JqGridUserRowConfig
from library.pandas_jqgrid import DataFrameJqGrid
from library.pandas_utils import nan_to_none
import pandas as pd
from seqauto.models import SeqAutoRun, SequencingSample, GoldCoverageSummary, \
    GoldReference, EnrichmentKit
from seqauto.seqauto_stats import get_sample_enrichment_kits_df, group_enrichment_kits_df


class SeqAutoRunsGrid(JqGridUserRowConfig):
    model = SeqAutoRun
    caption = 'SeqAutoRuns'

    fields = ['id', 'status', 'task_id',
              'created', 'scan_start', 'create_models_start', 'scripts_and_jobs_start', 'finish_date', 'job_launch_script_filename', 'error_exception']
    colmodel_overrides = {'id': {'width': 120, 'formatter': 'viewSeqAutoRunsLink'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class SequencingSamplesGrid(JqGridUserRowConfig):
    model = SequencingSample
    caption = 'Sequencing Samples'
    fields = ["sample_sheet__sequencing_run__name", "sample_name", "enrichment_kit__name"]
    colmodel_overrides = {'sample_sheet__sequencing_run__name': {'label': 'SequencingRun'},
                          'enrichment_kit__name': {'label': 'Enrichment Kit'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.get_current()
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("sequencing_samples_grid",
                                  kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'sample_sheet__sequencing_run__name',
                                  'sortorder': 'asc',
                                  'grid_export_url': grid_export_url})


class SequencingSamplesHistoricalGrid(DataFrameJqGrid):
    TIME_FRAMES = {"year": ("Sequencing Samples by year", 'Year'),
                   "month_offset": ("Sequencing Samples by month", 'YYMM'),
                   "all": ("All Sequencing Samples", 'Time Period')}

    def __init__(self, user, time_frame):
        super().__init__()
        if time_frame not in SequencingSamplesHistoricalGrid.TIME_FRAMES:
            msg = f"time_frame '{time_frame}' invalid! Must be one of: {','.join(SequencingSamplesHistoricalGrid.TIME_FRAMES)}"
            raise ValueError(msg)

        caption, column_label = SequencingSamplesHistoricalGrid.TIME_FRAMES[time_frame]
        self._overrides = {'ID': {'label': column_label}}
        self.caption = caption
        self.time_frame = time_frame
        grid_export_url = reverse("sequencing_samples_historical_grid",
                                  kwargs={"time_frame": time_frame,
                                          "op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'grid_export_url': grid_export_url})

    def get_dataframe(self):
        sample_enrichment_kits_df = get_sample_enrichment_kits_df()
        if self.time_frame == "all":
            start = sample_enrichment_kits_df["year_month"].min()
            end = sample_enrichment_kits_df["year_month"].max()
            enrichment_kit_counts = sample_enrichment_kits_df["enrichment_kit__name"].value_counts()
            enrichment_kit_counts.name = f"{start}-{end}"
            df = pd.DataFrame(enrichment_kit_counts).T
        else:
            enrichment_kits_over_time, enrichment_kit_labels = group_enrichment_kits_df(sample_enrichment_kits_df, self.time_frame)
            enrichment_kit_counts_dict = {k: v for (k, v) in enrichment_kits_over_time}
            df = pd.DataFrame.from_records(enrichment_kit_counts_dict, index=enrichment_kit_labels)

        return df


class GoldCoverageSummaryGrid(JqGridUserRowConfig):
    model = GoldCoverageSummary
    caption = 'GoldCoverageSummaries'
    fields = ["original_gene_symbol", "original_transcript_id", "gene_symbol", "transcript",
              "mean", "standard_error", "min_mean",
              "depth_20x_5th_percentile", "depth_10x_5th_percentile",
              "depth_mean_5th_percentile", "depth_mean_95th_percentile"]
    colmodel_overrides = {"gene_symbol": {"label": "Matched Symbol ", 'formatter': 'viewGoldSummaryGeneLink'},
                          "standard_error": {"server_side_formatter": lambda row, field: nan_to_none(row[field])}}

    def __init__(self, user, pk):
        super().__init__(user)
        gold_reference = get_object_or_404(GoldReference, pk=pk)
        queryset = self.model.objects.all()
        queryset = queryset.filter(gold_reference=gold_reference)

        self.queryset = queryset.values(*self.get_field_names())


class EnrichmentKitGeneCoverageGrid(JqGridUserRowConfig):
    SEQUENCING_SAMPLE_PATH = "gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample"
    SEQUENCING_RUN_PATH = SEQUENCING_SAMPLE_PATH + "__sample_sheet__sequencing_run"
    GOLD_PATH = SEQUENCING_RUN_PATH + "__gold_standard"
    SAMPLE_NAME_PATH = SEQUENCING_SAMPLE_PATH + "__sample_name"

    model = GeneCoverageCanonicalTranscript
    caption = 'Enrichment Kit Gene Coverage'
    fields = [SEQUENCING_RUN_PATH, GOLD_PATH, SAMPLE_NAME_PATH,
              "min", "mean", "std_dev", "percent_0x", "percent_10x", "percent_20x", "sensitivity"]
    number_format = {'formatter': 'number', 'width': 80}
    colmodel_overrides = {SEQUENCING_RUN_PATH: {"label": "Sequencing Run", "formatter": "viewSequencingRunLink"},
                          GOLD_PATH: {"label": "Gold"},
                          SAMPLE_NAME_PATH: {'label': 'Sample'},
                          'min': {'width': 40},
                          'mean': number_format,
                          'std_dev': number_format,
                          'percent_0x': number_format,
                          'percent_10x': number_format,
                          'percent_20x': number_format,
                          'sensitivity': number_format}

    def __init__(self, user, enrichment_kit_id, gene_symbol, extra_filters=None):
        super().__init__(user)
        enrichment_kit = get_object_or_404(EnrichmentKit, pk=enrichment_kit_id)
        if extra_filters is None:
            extra_filters = {}
        queryset = GeneCoverageCanonicalTranscript.filter_for_kit_and_gene_symbol(enrichment_kit, gene_symbol)
        if extra_filters.get("gold_only"):
            queryset = queryset.filter(**{self.GOLD_PATH: True})

        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': ",".join([self.SEQUENCING_SAMPLE_PATH, self.SAMPLE_NAME_PATH]),
                                  'sortorder': 'desc'})
