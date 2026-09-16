from django.db.models import QuerySet
from django.http import HttpRequest

from seqauto.models import FastQC, Flagstats, IlluminaFlowcellQC, QCExecSummary
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class IlluminaFlowcellQCColumns(DatatableConfig[IlluminaFlowcellQC]):
    grid_name = "IlluminaFlowcellQC"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="sample_sheet__sequencing_run__name", label="SequencingRun", orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key="sample_sheet__sequencing_run__sequencer__sequencer_model__model", label="Model",
                       orderable=True),
            RichColumn(key="sample_sheet__sequencing_run__sequencer", label="Sequencer", orderable=True),
            RichColumn(key="mean_cluster_density", label="Mean Cluster Density", orderable=True, css_class="num"),
            RichColumn(key="mean_pf_cluster_density", label="Mean PF Cluster Density", orderable=True, css_class="num"),
            RichColumn(key="total_clusters", label="Total Clusters", orderable=True, css_class="num"),
            RichColumn(key="total_pf_clusters", label="Total PF Clusters", orderable=True, css_class="num"),
            RichColumn(key="percentage_of_clusters_pf", label="% Clusters PF", orderable=True, css_class="num"),
            RichColumn(key="aligned_to_phix", label="Aligned To PhiX", orderable=True, css_class="num"),
        ]

    def get_initial_queryset(self) -> QuerySet[IlluminaFlowcellQC]:
        return IlluminaFlowcellQC.objects.all()


class FastQCColumns(DatatableConfig[FastQC]):
    grid_name = "FastQC"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="fastq__sequencing_sample__sample_sheet__sequencing_run__name", label="SequencingRun",
                       orderable=True, default_sort=SortOrder.DESC),
            RichColumn(key="fastq__name", label="FastQ", orderable=True),
            RichColumn(key="fastq__read", label="Read", orderable=True),
            RichColumn(key="total_sequences", label="Total Sequences", orderable=True, css_class="num"),
            RichColumn(key="filtered_sequences", label="Filtered Sequences", orderable=True, css_class="num"),
            RichColumn(key="gc", label="GC", orderable=True, css_class="num"),
        ]

    def get_initial_queryset(self) -> QuerySet[FastQC]:
        return FastQC.objects.all()


class FlagstatsColumns(DatatableConfig[Flagstats]):
    grid_name = "Flagstats"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="bam_file__sequencing_sample__sample_sheet__sequencing_run__name", label="SequencingRun",
                       orderable=True, default_sort=SortOrder.DESC),
            RichColumn(key="bam_file__name", label="BAM File", orderable=True),
            RichColumn(key="total", label="Total", orderable=True, css_class="num"),
            RichColumn(key="read1", label="Read1", orderable=True, css_class="num"),
            RichColumn(key="read2", label="Read2", orderable=True, css_class="num"),
            RichColumn(key="mapped", label="Mapped", orderable=True, css_class="num"),
            RichColumn(key="properly_paired", label="Properly Paired", orderable=True, css_class="num"),
        ]

    def get_initial_queryset(self) -> QuerySet[Flagstats]:
        return Flagstats.objects.all()


class QCExecSummaryColumns(DatatableConfig[QCExecSummary]):
    grid_name = "QCExecSummary"
    server_csv_download = True
    scroll_x = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="qc__bam_file__sequencing_sample__sample_sheet__sequencing_run__name",
                       label="SequencingRun", orderable=True, default_sort=SortOrder.DESC),
            RichColumn(key="qc__bam_file__sequencing_sample__sample_name", label="SampleName", orderable=True),
            RichColumn(key="percent_500x_goi", label="% 500x GOI", orderable=True, css_class="num"),
            RichColumn(key="percent_250x_goi", label="% 250x GOI", orderable=True, css_class="num"),
            RichColumn(key="percent_20x_goi", label="% 20x GOI", orderable=True, css_class="num"),
            RichColumn(key="percent_10x_goi", label="% 10x GOI", orderable=True, css_class="num"),
            RichColumn(key="mean_coverage_across_genes", label="Mean Coverage Across Genes", orderable=True,
                       css_class="num"),
            RichColumn(key="mean_coverage_across_kit", label="Mean Coverage Across Kit", orderable=True,
                       css_class="num"),
            RichColumn(key="uniformity_of_coverage", label="Uniformity Of Coverage", orderable=True, css_class="num"),
            RichColumn(key="percent_read_enrichment", label="% Read Enrichment", orderable=True, css_class="num"),
            RichColumn(key="median_insert", label="Median Insert", orderable=True, css_class="num"),
            RichColumn(key="ts_to_tv_ratio", label="Ts/Tv Ratio", orderable=True, css_class="num"),
            RichColumn(key="number_snps", label="Number SNPs", orderable=True, css_class="num"),
            RichColumn(key="snp_dbsnp_percent", label="SNP dbSNP %", orderable=True, css_class="num"),
            RichColumn(key="number_indels", label="Number Indels", orderable=True, css_class="num"),
            RichColumn(key="indels_dbsnp_percent", label="Indels dbSNP %", orderable=True, css_class="num"),
        ]

    def get_initial_queryset(self) -> QuerySet[QCExecSummary]:
        return QCExecSummary.objects.all()
