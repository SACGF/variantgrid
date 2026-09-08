import json
import logging
from collections import defaultdict

from django.db.models.aggregates import Count
from django.shortcuts import get_object_or_404, redirect, render

from library.utils import full_class_name
from seqauto import forms
from seqauto.graphs.index_metrics_qc_graph import IndexMetricsQCGraph
from seqauto.graphs.qc_exec_summary_graph import QCExecSummaryGraph
from seqauto.graphs.sequencing_run_qc_graph import SequencingRunQCGraph
from seqauto.models import (
    QC,
    EnrichmentKit,
    IlluminaFlowcellQC,
    QCColumn,
    QCExecSummary,
    QCType,
    SequencingRun,
    SequencingSample,
)
from seqauto.models.models_enums import QCCompareType
from seqauto.qc.sequencing_run_utils import (
    get_qc_exec_summary_data,
    get_sequencing_run_columns,
    get_sequencing_run_data,
)
from seqauto.seqauto_stats import get_sample_enrichment_kits_df
from snpdb.graphs import graphcache


def sequencing_stats(request):

    def count_values(klass, column):
        qs = klass.objects.values_list(column).annotate(Count("pk"))
        return list(qs)

    sr_sequencer = count_values(SequencingRun, "sequencer")
    sequencing_run_info = {
        'sequencer': sr_sequencer,
        'sequencer_model': count_values(SequencingRun, "sequencer__sequencer_model"),
        'enrichment_kit': count_values(SequencingRun, "enrichment_kit__name"),
    }

    ss_sequencer = count_values(SequencingSample, "sample_sheet__sequencing_run__sequencer")
    sequencing_sample_info = {
        'sequencer': ss_sequencer,
        'sequencer_model':  count_values(SequencingSample, "sample_sheet__sequencing_run__sequencer__sequencer_model"),
        'enrichment_kit': count_values(SequencingSample, "enrichment_kit__name"),
    }

    sample_enrichment_kits_df = get_sample_enrichment_kits_df()
    context = {'num_sequencing_runs': sum(i[1] for i in sr_sequencer),
               'num_samples': sum(i[1] for i in ss_sequencer),
               'sequencing_run_info': sequencing_run_info,
               'sequencing_sample_info': sequencing_sample_info,
               'sample_enrichment_kits_df': sample_enrichment_kits_df}
    return render(request, 'seqauto/sequencing_stats.html', context)


def sequencing_stats_data(request):
    context = {"has_sequencing_samples": SequencingSample.objects.all().exists()}
    return render(request, 'seqauto/sequencing_stats_data.html', context)


def qc_data(request):
    return render(request, 'seqauto/qc_data.html')


def qc_graphs(request):
    form = forms.QCColumnForm()
    qc_type_totals = dict(QCType.objects.all().values_list("name", "total_field"))

    context = {
        'form': form,
        'qc_type_totals': qc_type_totals,
    }
    return render(request, 'seqauto/qc_graphs.html', context)


def sequencing_run_qc_graph(request, sequencing_run_id, qc_compare_type):
    _ = QCCompareType(qc_compare_type)  # Check valid
    graph_class_name = full_class_name(SequencingRunQCGraph)
    cached_graph = graphcache.async_graph(graph_class_name, sequencing_run_id, qc_compare_type)
    return redirect(cached_graph)


def sequencing_run_qc_json_graph(request, sequencing_run_id, qc_compare_type):
    qc_compare_types_dict = dict(QCCompareType.choices)
    if qc_compare_type not in qc_compare_types_dict:
        msg = f"Unknown QCCompareType '{qc_compare_type}'"
        raise ValueError(msg)

    sequencing_run = SequencingRun.objects.get(pk=sequencing_run_id)
    sequencing_run_data = get_sequencing_run_data(sequencing_run, qc_compare_type, include_passed_sequencing_run=True)
    ss_path = IlluminaFlowcellQC.get_sequencing_run_path()
    sr_columns = get_sequencing_run_columns(ss_path, ['name', 'gold_standard'])
    (sequencing_run_column, gold_column) = tuple(sr_columns)

    context = {
        "container_name": "run-stats",
        "current_label": sequencing_run_id,
        "qc_data": sequencing_run_data,
        "label_column": sequencing_run_column,
        "gold_column": gold_column
    }
    return render(request, 'seqauto/json_graphs/qc_json_graph.html', context)


def _qc_column_box_data(qc_column: QCColumn, use_percent: bool) -> dict[str, list]:
    """ QC values per enrichment kit (full name), optionally as a percentage of the QC type's total field """
    sequencing_sample_path = 'bam_file__sequencing_sample'
    enrichment_kit_path = sequencing_sample_path + '__enrichment_kit'

    def get_field(f):
        return qc_column.qc_type.qc_object_path + "__" + f

    path = get_field(qc_column.field)
    qs = QC.objects.filter(**{path + "__isnull": False})
    qs = qs.order_by(enrichment_kit_path)  # want dict eventually sorted by kit

    total_field = None
    args = [enrichment_kit_path + "__name", enrichment_kit_path + "__version", path]
    if use_percent:
        total_field = qc_column.qc_type.total_field
        if total_field is None:
            msg = f"Asked for percentage for {qc_column} ({qc_column.qc_type}) with no total field!"
            raise ValueError(msg)
        total_path = get_field(total_field)
        args.append(total_path)

    kit_values = defaultdict(list)
    for values in qs.values_list(*args):
        enrichment_kit_name = values[0]
        enrichment_kit_version = values[1]
        name = EnrichmentKit.get_full_name(enrichment_kit_name, enrichment_kit_version)
        if use_percent:
            (val, total) = values[2:]
            if total:
                val = 100.0 * val / total
            elif val:
                msg = f"Val was '{val}' with total field {total_field} of {total}"
                raise ValueError(msg)
        else:
            val = values[2]
        kit_values[name].append(val)
    return dict(kit_values)


def qc_column_graph(request, qc_column_id, use_percent):
    qc_column = get_object_or_404(QCColumn, pk=qc_column_id)
    logging.info("Using %s", qc_column)
    use_percent = json.loads(use_percent)  # Boolean

    graph_title = str(qc_column)
    if use_percent:
        graph_title += " %"

    context = {
        "box_data": _qc_column_box_data(qc_column, use_percent),
        "graph_title": graph_title,
    }
    return render(request, 'seqauto/json_graphs/qc_column_graph.html', context)


def index_metrics_qc_graph(request, illumina_qc_id):
    graph_class_name = full_class_name(IndexMetricsQCGraph)
    cached_graph = graphcache.async_graph(graph_class_name, illumina_qc_id)
    return redirect(cached_graph)


def qc_exec_summary_graph(request, qc_exec_summary_id, qc_compare_type):
    graph_class_name = full_class_name(QCExecSummaryGraph)
    cached_graph = graphcache.async_graph(graph_class_name, qc_exec_summary_id, qc_compare_type)
    return redirect(cached_graph)


def qc_exec_summary_json_graph(request, qc_exec_summary_id, qc_compare_type):
    _ = QCCompareType(qc_compare_type)  # Check valid

    def get_label(sequencing_run_name, sample_name):
        return f"{sequencing_run_name}__newline__{sample_name}"

    qc_exec_summary = QCExecSummary.objects.get(pk=qc_exec_summary_id)
    qc_exec_summary_data = get_qc_exec_summary_data(qc_exec_summary.sequencing_run, qc_compare_type, qc_exec_summary, include_passed_qc_exec_summary=True)
    ss_path = QCExecSummary.get_sequencing_run_path()
    sr_columns = get_sequencing_run_columns(ss_path, ['name', 'gold_standard'])
    (sequencing_run_column, gold_column) = tuple(sr_columns)

    # Create a new label based on sequencing_run + sample
    current_label = get_label(qc_exec_summary.sequencing_run.name, qc_exec_summary.sample_name)
    sequencing_run_names = qc_exec_summary_data[sequencing_run_column]
    sample_names = qc_exec_summary_data["qc__bam_file__sequencing_sample__sample_name"]
    labels = [get_label(sr, ss) for sr, ss in zip(sequencing_run_names, sample_names)]
    qc_exec_summary_data["label"] = labels

    context = {
        "container_name": "exec-summary",
        "current_label": current_label,
        "qc_data": qc_exec_summary_data,
        "label_column": "label",
        "gold_column": gold_column
    }
    return render(request, 'seqauto/json_graphs/qc_json_graph.html', context)
