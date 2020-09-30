from django.urls.base import reverse

from library.django_utils.jqgrid_view import JQGridViewOp
from library.jqgrid_user_row_config import JqGridUserRowConfig
from seqauto.models import IlluminaFlowcellQC, FastQC, Flagstats, QCExecSummary


class IlluminaFlowcellQCGrid(JqGridUserRowConfig):
    model = IlluminaFlowcellQC
    caption = 'IlluminaFlowcellQC'
    fields = ["sample_sheet__sequencing_run__name", "sample_sheet__sequencing_run__sequencer__sequencer_model__model", "sample_sheet__sequencing_run__sequencer",
              "mean_cluster_density", "mean_pf_cluster_density", "total_clusters", "total_pf_clusters", "percentage_of_clusters_pf", "aligned_to_phix"]
    colmodel_overrides = {'sample_sheet__sequencing_run__name': {'label': 'SequencingRun'},
                          'sample_sheet__sequencing_run__sequencer__sequencer_model__model': {'label': 'Model'},
                          'sample_sheet__sequencing_run__sequencer': {'label': 'Sequencer'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.filter(data_state='C')
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("illumina_flowcell_qc_grid", kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'sample_sheet__sequencing_run__name',
                                  'sortorder': 'desc',
                                  'grid_export_url': grid_export_url})


class FastQCGrid(JqGridUserRowConfig):
    model = FastQC
    caption = 'FastQC'
    fields = ["fastq__sequencing_sample__sample_sheet__sequencing_run__name", "fastq__name", "fastq__read", "total_sequences", "filtered_sequences", "gc"]
    colmodel_overrides = {'fastq__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'SequencingRun'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.filter(data_state='C')
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("fastqc_grid", kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'fastq__sequencing_sample__sample_sheet__sequencing_run__name',
                                  'sortorder': 'desc',
                                  'grid_export_url': grid_export_url})


class FlagstatsGrid(JqGridUserRowConfig):
    model = Flagstats
    caption = 'Flagstats'
    fields = ["bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name",
              "bam_file__name", "total", "read1", "read2", "mapped", "properly_paired"]
    colmodel_overrides = {'bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'SequencingRun'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.filter(data_state='C')
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("flagstats_grid", kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name',
                                  'sortorder': 'desc',
                                  'grid_export_url': grid_export_url})


class QCExecSummaryGrid(JqGridUserRowConfig):
    model = QCExecSummary
    caption = 'QCExecSummary'
    fields = ["qc__bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name",
              "qc__bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_name",
              "percent_500x", "percent_250x", "percent_20x", "percent_10x", "mean_coverage_across_genes", "mean_coverage_across_kit",
              "uniformity_of_coverage", "percent_read_enrichment", "duplicated_alignable_reads", "median_insert", "ts_to_tv_ratio",
              "number_snps", "snp_dbsnp_percent", "number_indels", "indels_dbsnp_percent"]
    colmodel_overrides = {'qc__bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name': {'label': 'SequencingRun'},
                          "qc__bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_name": {'label': 'SampleName'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.objects.filter(data_state='C')
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("qc_exec_summary_grid", kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'qc__bam_file__unaligned_reads__fastq_r1__sequencing_sample__sample_sheet__sequencing_run__name',
                                  'sortorder': 'desc',
                                  'grid_export_url': grid_export_url})
