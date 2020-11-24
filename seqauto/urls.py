from rest_framework.urlpatterns import format_suffix_patterns

from library.django_utils.jqgrid_view import JQGridView
from seqauto import views, views_autocomplete, views_rest
from seqauto.grids.qc_data_grids import IlluminaFlowcellQCGrid, FastQCGrid, FlagstatsGrid, \
    QCExecSummaryGrid
from seqauto.grids.seqauto_grids import SeqAutoRunsGrid, EnrichmentKitGeneCoverageGrid, \
    GoldCoverageSummaryGrid, SequencingSamplesGrid, SequencingSamplesHistoricalGrid
from seqauto.grids.sequencing_data_grids import ExperimentsListGrid, SequencingRunListGrid, \
    UnalignedReadsListGrid, BamFileListGrid, VCFFileListGrid, QCFileListGrid, \
    EnrichmentKitListGrid
from seqauto.grids.sequencing_software_versions_grids import LibraryGrid, SequencerGrid, \
    AssayGrid, AlignerGrid, VariantCallerGrid, VariantCallingPipelineGrid
from seqauto.views import SequencerUpdate, LibraryUpdate, AssayUpdate, VariantCallerUpdate, \
    AlignerUpdate, VariantCallingPipelineUpdate
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('', views.sequencing_data, name='sequencing_data'),
    perm_path('seqauto_runs', views.seqauto_runs, name='seqauto_runs'),
    perm_path('experiments', views.experiments, name='experiments'),
    perm_path('sequencing_runs', views.sequencing_runs, name='sequencing_runs'),
    perm_path('unaligned_reads', views.unaligned_reads, name='unaligned_reads'),
    perm_path('bam_files', views.bam_files, name='bam_files'),
    perm_path('vcf_files', views.vcf_files, name='vcf_files'),
    perm_path('qcs', views.qcs, name='qcs'),
    perm_path('view_experiment/<experiment_id>', views.view_experiment, name='view_experiment'),
    perm_path('software_pipeline', views.software_pipeline, name='software_pipeline'),

    perm_path('enrichment_kits_list', views.enrichment_kits_list, name='enrichment_kits_list'),
    perm_path('enrichment_kit/<int:pk>', views.view_enrichment_kit, name='view_enrichment_kit'),

    perm_path('sequencing_stats', views.sequencing_stats, name='sequencing_stats'),
    perm_path('sequencing_stats/data', views.sequencing_stats_data, name='sequencing_stats_data'),
    perm_path('qc_data', views.qc_data, name='qc_data'),
    perm_path('sequencing_historical_graphs', views.sequencing_historical_graphs, name='sequencing_historical_graphs'),
    perm_path('graphs/qc_column_historical_graph/<int:qc_column_id>/<graph_type>/<enrichment_kit_separation>/<enrichment_kit_id>/<use_percent>', views.qc_column_historical_graph, name='qc_column_historical_graph'),
    perm_path('graphs/sequencing_run_qc_graph/<sequencing_run_id>/<qc_compare_type>', views.sequencing_run_qc_graph, name='sequencing_run_qc_graph'),
    perm_path('graphs/sequencing_run_qc_json_graph/<sequencing_run_id>/<qc_compare_type>', views.sequencing_run_qc_json_graph, name='sequencing_run_qc_json_graph'),

    perm_path('graphs/index_metrics_qc_graph/<illumina_qc_id>', views.index_metrics_qc_graph, name='index_metrics_qc_graph'),
    perm_path('graphs/qc_exec_summary_graph/<qc_exec_summary_id>/<qc_compare_type>', views.qc_exec_summary_graph, name='qc_exec_summary_graph'),
    perm_path('graphs/qc_exec_summary_json_graph/<qc_exec_summary_id>/<qc_compare_type>', views.qc_exec_summary_json_graph, name='qc_exec_summary_json_graph'),

    perm_path('view_seqauto_run/<int:seqauto_run_id>', views.view_seqauto_run, name='view_seqauto_run'),
    perm_path('view_sequencing_run/<sequencing_run_id>/tab/<int:tab_id>', views.view_sequencing_run, name='view_sequencing_run_tab'),
    perm_path('view_sequencing_run_stats_tab/<sequencing_run_id>', views.view_sequencing_run_stats_tab, name='view_sequencing_run_stats_tab'),
    perm_path('view_sequencing_run/<sequencing_run_id>', views.view_sequencing_run, name='view_sequencing_run'),
    perm_path('sequencing_run/reload_experiment_name/<sequencing_run_id>', views.reload_experiment_name, name='reload_experiment_name'),
    perm_path('sequencing_run/delete/<sequencing_run_id>', views.delete_sequencing_run, name='delete_sequencing_run'),
    perm_path('sequencing_run/assign_data_to_current_sample_sheet/<sequencing_run_id>', views.assign_data_to_current_sample_sheet, name='assign_data_to_current_sample_sheet'),

    perm_path('view_unaligned_reads/<int:unaligned_reads_id>', views.view_unaligned_reads, name='view_unaligned_reads'),
    perm_path('view_bam/<int:bam_file_id>', views.view_bam_file, name='view_bam_file'),
    perm_path('view_combo_vcf_file/<int:combo_vcf_file_id>', views.view_combo_vcf_file, name='view_combo_vcf_file'),
    perm_path('view_vcf/<int:vcf_file_id>', views.view_vcf_file, name='view_vcf_file'),
    perm_path('view_qc/<int:qc_id>', views.view_qc, name='view_qc'),
    perm_path('view_qc/view_qc_exec_summary_tab/<int:qc_id>', views.view_qc_exec_summary_tab, name='view_qc_exec_summary_tab'),
    perm_path('view_qc/view_qc_gene_list_tab/<int:qc_id>', views.view_qc_gene_list_tab, name='view_qc_gene_list_tab'),
    perm_path('view_qc/view_gene_coverage_collection_tab/<int:qc_id>', views.view_qc_gene_coverage_collection_tab, name='view_qc_gene_coverage_collection_tab'),
    perm_path('view_sample_qc_tab/<int:sample_id>', views.view_sample_qc_tab, name='view_sample_qc_tab'),
    perm_path('view_enrichment_kit_gene_coverage/<int:enrichment_kit_id>/<gene_symbol>', views.view_enrichment_kit_gene_coverage, name='view_enrichment_kit_gene_coverage'),
    perm_path('view_gold_coverage_summary/<int:pk>', views.view_gold_coverage_summary, name='view_gold_coverage_summary'),

    # Grids
    perm_path('seqauto_runs/grid/<slug:op>/', JQGridView.as_view(grid=SeqAutoRunsGrid), name='seqauto_runs_grid'),
    perm_path('experiments/grid/<slug:op>/', JQGridView.as_view(grid=ExperimentsListGrid, delete_row=True), name='experiments_grid'),
    perm_path('sequencing_run/grid/<slug:op>/', JQGridView.as_view(grid=SequencingRunListGrid), name='sequencing_run_grid'),
    perm_path('unaligned_reads/grid/<slug:op>/', JQGridView.as_view(grid=UnalignedReadsListGrid), name='unaligned_reads_grid'),
    perm_path('bam_file/grid/<slug:op>/', JQGridView.as_view(grid=BamFileListGrid), name='bam_file_grid'),
    perm_path('vcf_file/grid/<slug:op>/', JQGridView.as_view(grid=VCFFileListGrid), name='vcf_file_grid'),
    perm_path('qc/grid/<slug:op>/', JQGridView.as_view(grid=QCFileListGrid), name='qc_grid'),
    perm_path('enrichment_kit/grid/<slug:op>/', JQGridView.as_view(grid=EnrichmentKitListGrid), name='enrichment_kit_list_grid'),
    perm_path('enrichment_kit/gene/grid/<int:enrichment_kit_id>/<genome_build_name>/<gene_symbol>/<slug:op>/', JQGridView.as_view(grid=EnrichmentKitGeneCoverageGrid), name='enrichment_kit_gene_coverage_grid'),
    perm_path('gold_coverage_summary/grid/<pk>/<slug:op>/', JQGridView.as_view(grid=GoldCoverageSummaryGrid), name='gold_coverage_summary_grid'),
    perm_path('sequencing_stats/sequencing_samples/grid/<slug:op>/', JQGridView.as_view(grid=SequencingSamplesGrid, csv_download=True), name='sequencing_samples_grid'),
    perm_path('sequencing_stats/sequencing_samples_historical/grid/<slug:time_frame>/<slug:op>/', JQGridView.as_view(grid=SequencingSamplesHistoricalGrid, csv_download=True), name='sequencing_samples_historical_grid'),
    # QC Data
    perm_path('illumina_flowcell_qc/grid/<slug:op>/', JQGridView.as_view(grid=IlluminaFlowcellQCGrid, csv_download=True), name='illumina_flowcell_qc_grid'),
    perm_path('fastqc/grid/<slug:op>/', JQGridView.as_view(grid=FastQCGrid, csv_download=True), name='fastqc_grid'),
    perm_path('flagstats/grid/<slug:op>/', JQGridView.as_view(grid=FlagstatsGrid, csv_download=True), name='flagstats_grid'),
    perm_path('qc_exec_summary/grid/<slug:op>/', JQGridView.as_view(grid=QCExecSummaryGrid, csv_download=True), name='qc_exec_summary_grid'),
    # Software/settings
    perm_path('sequencing_software_versions/library/grid/<slug:op>/', JQGridView.as_view(grid=LibraryGrid), name='library_grid'),
    perm_path('sequencing_software_versions/sequencer/grid/<slug:op>/', JQGridView.as_view(grid=SequencerGrid), name='sequencer_grid'),
    perm_path('sequencing_software_versions/assay/grid/<slug:op>/', JQGridView.as_view(grid=AssayGrid), name='assay_grid'),
    perm_path('sequencing_software_versions/aligner/grid/<slug:op>/', JQGridView.as_view(grid=AlignerGrid), name='aligner_grid'),
    perm_path('sequencing_software_versions/variant_caller/grid/<slug:op>/', JQGridView.as_view(grid=VariantCallerGrid), name='variant_caller_grid'),
    perm_path('sequencing_software_versions/variant_calling_pipeline/grid/<slug:op>/', JQGridView.as_view(grid=VariantCallingPipelineGrid), name='variant_calling_pipeline_grid'),

    perm_path('sequencing_software_versions', views.sequencing_software_versions, name='sequencing_software_versions'),
    perm_path('view_sequencer/<pk>', SequencerUpdate.as_view(), name='view_sequencer'),
    perm_path('view_library/<pk>', LibraryUpdate.as_view(), name='view_library'),
    perm_path('view_assay/<pk>', AssayUpdate.as_view(), name='view_assay'),
    perm_path('view_variant_caller/<pk>', VariantCallerUpdate.as_view(), name='view_variant_caller'),
    perm_path('view_aligner/<pk>', AlignerUpdate.as_view(), name='view_aligner'),
    perm_path('view_variant_calling_pipeline/<pk>', VariantCallingPipelineUpdate.as_view(), name='view_variant_calling_pipeline'),

    # Autocompletes
    perm_path('autocomplete/QCColumn/', views_autocomplete.QCColumnAutocompleteView.as_view(), name='qc_column_autocomplete'),
    perm_path('autocomplete/EnrichmentKit/', views_autocomplete.EnrichmentKitAutocompleteView.as_view(), name='enrichment_kit_autocomplete'),
]

#router = routers.DefaultRouter()
rest_urlpatterns = [
    perm_path('api/view_enrichment_kit_summary/<int:pk>', views_rest.EnrichmentKitSummaryView.as_view(), name='api_view_enrichment_kit_summary'),
    perm_path('api/view_enrichment_kit/<int:pk>', views_rest.EnrichmentKitView.as_view(), name='api_view_enrichment_kit'),
    perm_path('api/enrichment_kit_gene_coverage/<int:enrichment_kit_id>/<gene_symbol>', views_rest.EnrichmentKitGeneCoverageView.as_view(), name='api_enrichment_kit_gene_coverage'),
    perm_path('api/enrichment_kit_gene_gold_coverage/<int:enrichment_kit_id>/<gene_symbol>', views_rest.EnrichmentKitGeneGoldCoverageView.as_view(), name='api_enrichment_kit_gene_gold_coverage'),
    perm_path('api/enrichment_kit_gene_gold_coverage_summary/<int:enrichment_kit_id>/<gene_symbol>', views_rest.GoldCoverageSummaryView.as_view(), name='api_enrichment_kit_gene_gold_coverage_summary'),
    perm_path('api/enrichment_kit_gene_gold_coverage_summary/batch/<int:enrichment_kit_id>', views_rest.BatchGoldCoverageSummaryView.as_view(), name='api_batch_enrichment_kit_gene_gold_coverage_summary'),
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
