"""
Created on 30/05/2016

Place to put common data gathering utils for SequencingRunQC and QCExecsummary
(ie used by both Matplotlib and JSON QC graphs)

"""

from collections import defaultdict

from seqauto.models import SequencingRun, IlluminaFlowcellQC, ReadQ30, QCExecSummary, QCType
from seqauto.models.models_enums import QCCompareType


ILLUMINA_FLOWCELL_QC_COLUMNS = ["mean_cluster_density", "mean_pf_cluster_density", "total_clusters",
                                "total_pf_clusters", "percentage_of_clusters_pf", "aligned_to_phix"]
PAIRED_END_READS = ('R1', 'R2')


def get_q30_col_name(read):
    return "%s %% >= Q30" % read


READ_COLUMNS = [get_q30_col_name(read) for read in PAIRED_END_READS]
SEQUENCING_RUN_QC_COLUMNS = ILLUMINA_FLOWCELL_QC_COLUMNS + READ_COLUMNS


def get_sequencing_runs_for_compare_type(sequencing_run, qc_compare_type):
    all_runs_qs = SequencingRun.objects.filter(bad=False)

    if qc_compare_type == QCCompareType.ALL_RUNS:
        qs = all_runs_qs
    elif qc_compare_type == QCCompareType.SEQUENCER:
        qs = all_runs_qs.filter(sequencer=sequencing_run.sequencer)
    elif qc_compare_type in QCCompareType.ENRICHMENT_KIT_TYPES:
        enrichment_kit = sequencing_run.enrichment_kit
        qs = all_runs_qs.filter(seqautorecord__samplesheet__sequencingsample__enrichment_kit=enrichment_kit).distinct()
        if qc_compare_type == QCCompareType.GOLD_ENRICHMENT_KIT_RUNS:
            qs = qs.filter(gold_standard=True)
    elif qc_compare_type == QCCompareType.SEQUENCING_RUN:
        qs = all_runs_qs.filter(pk=sequencing_run.pk)
    else:
        msg = f"SequencingRunQCGraph: Unknown qc_compare_type '{qc_compare_type}'"
        raise ValueError(msg)

    return qs


def get_sequencing_runs(sequencing_run, qc_compare_type, include_passed_sequencing_run):
    """ Used to get 'other' sequencing runs for QCCompareType
        but pass include_passed_sequencing_run to also include this """
    sequencing_runs_qs = get_sequencing_runs_for_compare_type(sequencing_run, qc_compare_type)
    if not include_passed_sequencing_run:
        sequencing_runs_qs = sequencing_runs_qs.exclude(pk=sequencing_run.pk)

    sequencing_runs_ids = set(sequencing_runs_qs.values_list("pk", flat=True))
    if include_passed_sequencing_run:
        # Need to work with sets not QS as get error "Cannot combine a unique query with a non-unique query"
        sequencing_runs_ids.add(sequencing_run.pk)

    return sequencing_runs_ids


def get_sequencing_run_columns(ss_path, columns):
    return [f"{ss_path}__{c}" for c in columns]


def get_sequencing_run_data(sequencing_run, qc_compare_type, include_passed_sequencing_run=False):
    sequencing_runs_ids = get_sequencing_runs(sequencing_run, qc_compare_type, include_passed_sequencing_run)
    run_data = defaultdict(list)
    flowcell_qc_qs = IlluminaFlowcellQC.objects.filter(data_state='C', sample_sheet__sequencing_run__in=sequencing_runs_ids)
    flowcell_qc_qs = flowcell_qc_qs.order_by("sample_sheet__sequencing_run__name")
    ss_path = IlluminaFlowcellQC.get_sequencing_run_path()
    values = ILLUMINA_FLOWCELL_QC_COLUMNS + get_sequencing_run_columns(ss_path, ['name', 'gold_standard'])
    non_null_kwargs = {"%s__isnull" % f: False for f in values}

    for data in flowcell_qc_qs.filter(**non_null_kwargs).values(*values):
        for k, v in data.items():
            run_data[k].append(v)

    read_q30_qs = ReadQ30.objects.filter(illumina_flowcell_qc__sample_sheet__sequencing_run__in=sequencing_runs_ids, read__in=PAIRED_END_READS)
    read_q30_qs = read_q30_qs.order_by("illumina_flowcell_qc__sample_sheet__sequencing_run__name")
    for read, percent in read_q30_qs.values_list("read", "percent"):
        run_data[get_q30_col_name(read)].append(percent)

    return run_data


def get_qc_exec_summary_data(sequencing_run, qc_compare_type, qc_exec_summary, include_passed_qc_exec_summary=False):
    sequencing_runs_ids = get_sequencing_runs(sequencing_run, qc_compare_type, True)
    ss_path = QCExecSummary.get_sequencing_run_path()
    run_data = defaultdict(list)

    kwargs = {"data_state": 'C',
              ss_path + "__in": sequencing_runs_ids}
    qc_exec_qs = QCExecSummary.objects.filter(**kwargs)
    if not include_passed_qc_exec_summary:
        qc_exec_qs = qc_exec_qs.exclude(pk=qc_exec_summary.pk)

    exec_summary_qc = QCType.objects.get(name="ExecSummaryQC")
    qc_exec_summary_columns = list(exec_summary_qc.qccolumn_set.all().values_list("field", flat=True))

    coverage_columns = list(qc_exec_summary.get_coverage_columns())
    sequencing_sample = "qc__bam_file__unaligned_reads__sequencing_sample__sample_name"
    sequencing_run_columns = get_sequencing_run_columns(ss_path, ['name', 'gold_standard'])
    values = ["pk", sequencing_sample] + qc_exec_summary_columns + coverage_columns + sequencing_run_columns
    non_null_kwargs = {"%s__isnull" % f: False for f in coverage_columns}

    for data in qc_exec_qs.filter(**non_null_kwargs).values(*values):
        for k, v in data.items():
            run_data[k].append(v)

    return run_data
