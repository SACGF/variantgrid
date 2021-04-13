import json
import os
from collections import defaultdict

from django.conf import settings
from django.contrib import messages
from django.core.exceptions import PermissionDenied
from django.db.models.aggregates import Count
from django.http.response import HttpResponseRedirect, JsonResponse, HttpResponse
from django.shortcuts import render, get_object_or_404
from django.urls.base import reverse
from django.utils.decorators import method_decorator
from django.views.decorators.http import require_POST
from django.views.generic.edit import UpdateView

from eventlog.models import create_event
from genes.models import CanonicalTranscriptCollection
from library.django_utils import get_model_fields, staff_only
from library.log_utils import log_traceback
from library.utils import full_class_name
from seqauto import forms
from seqauto.forms import SequencingRunForm, AllEnrichmentKitForm
from seqauto.graphs.index_metrics_qc_graph import IndexMetricsQCGraph
from seqauto.graphs.qc_column_boxplot_graph import QCColumnBoxPlotGraph
from seqauto.graphs.qc_column_line_graph import QCColumnLineGraph
from seqauto.graphs.qc_exec_summary_graph import QCExecSummaryGraph
from seqauto.graphs.sequencing_run_qc_graph import SequencingRunQCGraph
from seqauto.illumina.run_parameters import get_run_parameters
from seqauto.models import BamFile, SequencingRun, FastQC, Flagstats, UnalignedReads, QCType, VCFFile, QC, \
    Experiment, SequencingSample, SampleSheetCombinedVCFFile, QCExecSummary, IlluminaFlowcellQC, SeqAutoRun, \
    Library, Sequencer, Assay, Aligner, VariantCaller, VariantCallingPipeline, SoftwarePipelineNode, \
    GoldReference, GoldGeneCoverageCollection, EnrichmentKit, QCGeneCoverage
from seqauto.models.models_enums import QCGraphEnrichmentKitSeparationChoices, QCGraphType, \
    QCCompareType, SequencingFileType
from seqauto.qc.sequencing_run_utils import get_sequencing_run_data, get_qc_exec_summary_data, \
    get_sequencing_run_columns, SEQUENCING_RUN_QC_COLUMNS, QC_EXEC_SUMMARY_QC_COLUMNS
from seqauto.seqauto_stats import get_sample_enrichment_kits_df
from seqauto.sequencing_files.create_resource_models import assign_old_sample_sheet_data_to_current_sample_sheet
from seqauto.tasks.scan_run_jobs import scan_run_jobs
from snpdb.graphs import graphcache
from snpdb.models import Sample, UserSettings


def sequencing_data(request):
    context = {"last_success_datetime": SeqAutoRun.get_last_success_datetime()}
    return render(request, 'seqauto/sequencing_data.html', context)


def seqauto_runs(request):
    # Only allow button if settings allow, and previous run has finished
    enable_button = request.user.is_staff or settings.SEQAUTO_ALLOW_NON_STAFF_MANUAL_RUN

    if request.method == "POST":
        if not enable_button:
            msg = "Not allowed to manually kick off sequencing run"
            raise PermissionDenied(msg)

        only_process_file_types = []  # All
        only_launch_file_types = [SequencingFileType.ILLUMINA_FLOWCELL_QC]

        task = scan_run_jobs.si(only_process_file_types=only_process_file_types,  # @UndefinedVariable
                                only_launch_file_types=only_launch_file_types)
        task.apply_async()

        msg = 'Scanning disk for sequencing data...'
        messages.add_message(request, messages.INFO, msg, extra_tags='import-message')
        enable_button = False

    context = {"last_success_datetime": SeqAutoRun.get_last_success_datetime(),
               "enable_button": enable_button}
    return render(request, 'seqauto/seqauto_runs.html', context)


def experiments(request):
    return render(request, 'seqauto/experiments.html')


def sequencing_runs(request):
    context = {
        "enrichment_kit_form": AllEnrichmentKitForm(),
    }
    return render(request, 'seqauto/sequencing_runs.html', context)


def unaligned_reads(request):
    return render(request, 'seqauto/unaligned_reads.html')


def bam_files(request):
    return render(request, 'seqauto/bam_files.html')


def vcf_files(request):
    return render(request, 'seqauto/vcf_files.html')


def qcs(request):
    return render(request, 'seqauto/qc.html')


def view_experiment(request, experiment_id):
    experiment = get_object_or_404(Experiment, pk=experiment_id)
    context = {"experiment": experiment}
    return render(request, 'seqauto/view_experiment.html', context)


def software_pipeline(request):
    """ Extract sofware pipelines into Directed Acyclic Graphs """

    qs = SoftwarePipelineNode.objects.filter()  # @UndefinedVariable
    return render(request, 'seqauto/software_pipeline.html', {'dag_list': qs})


def view_seqauto_run(request, seqauto_run_id):

    seqauto_run = get_object_or_404(SeqAutoRun, pk=seqauto_run_id)
    qs = seqauto_run.jobscript_set.all()
    file_type_counts = qs.values("file_type").annotate(file_type_count=Count("file_type"))
    pbs_script_counts = []

    for file_type, file_type_count in file_type_counts.values_list("file_type", "file_type_count"):
        pbs_script_counts.append((SequencingFileType(file_type).label, file_type_count))

    if seqauto_run.error_exception:
        status = messages.ERROR
        messages.add_message(request, status, seqauto_run.error_exception, extra_tags='import-message')

    context = {"seqauto_run": seqauto_run,
               "pbs_script_counts": pbs_script_counts}
    return render(request, 'seqauto/view_seqauto_run.html', context)


def get_illumina_qc_and_show_stats_for_sample_sheet(sample_sheet):
    try:
        illumina_qc = sample_sheet.illuminaflowcellqc
        show_stats = illumina_qc.data_state == 'C'
    except:
        illumina_qc = None
        show_stats = False

    return illumina_qc, show_stats


def view_sequencing_run(request, sequencing_run_id, tab_id=0):
    sequencing_run = get_object_or_404(SequencingRun, pk=sequencing_run_id)

    sequencing_run_form = SequencingRunForm(request.POST or None, instance=sequencing_run)
    if request.method == "POST":
        if not request.user.is_staff:
            msg = "Only admin users can alter sequencing runs!"
            raise PermissionDenied(msg)

        if sequencing_run_form.is_valid():
            for f in sequencing_run_form.changed_data:
                val = sequencing_run_form.cleaned_data.get(f)
                message = f"{f} set to {val}"
                #SequencingRunModification.objects.create(sequencing_run=sequencing_run,
                #                                         user=request.user,
                #                                         message=message)

            sequencing_run = sequencing_run_form.save()

    if not sequencing_run.is_valid:  # Had errors
        sequencing_run.save()  # Try again now
        if sequencing_run.is_valid:
            message = "SequencingRun had errors but they appear to have been resolved. Setting is_valid=True"
            messages.add_message(request, messages.WARNING, message)

    sequencing_run.add_messages(request)

    run_vcfs = []
    try:
        for vcf_for_run in sequencing_run.vcffromsequencingrun_set.all():
            run_vcfs.append((vcf_for_run.variant_caller, vcf_for_run.vcf, vcf_for_run.vcf.can_view(request.user)))
    except:
        log_traceback()
    context = {
        "sequencing_run": sequencing_run,
        "sequencing_run_form": sequencing_run_form,
        'tab_id': tab_id,
        'run_vcfs': run_vcfs
    }

    try:  # May not have sample sheet and die
        sample_sheet = sequencing_run.get_current_sample_sheet()
        illumina_qc, show_stats = get_illumina_qc_and_show_stats_for_sample_sheet(sample_sheet)
        has_sequencing_sample_data = sample_sheet.sequencingsample_set.filter(sequencingsampledata__isnull=False).exists()
        context['sample_sheet'] = sample_sheet
        context['illumina_qc'] = illumina_qc
        context['show_stats'] = show_stats
        context['has_sequencing_sample_data'] = has_sequencing_sample_data
        context["sequencing_samples"] = sample_sheet.get_sorted_sequencing_samples()
        context['data_out_of_date_from_current_sample_sheet'] = sequencing_run.is_data_out_of_date_from_current_sample_sheet
    except:
        log_traceback()

    return render(request, 'seqauto/view_sequencing_run.html', context)


def view_sequencing_run_stats_tab(request, sequencing_run_id):
    sequencing_run = get_object_or_404(SequencingRun, pk=sequencing_run_id)

    compare_against = [QCCompareType.ALL_RUNS]
    if sequencing_run.sequencer.sequencingrun_set.exclude(pk=sequencing_run).exists():
        compare_against.append(QCCompareType.SEQUENCER)

    if sequencing_run.enrichment_kit:
        if sequencing_run.enrichment_kit.sequencingrun_set.exclude(pk=sequencing_run).exists():
            compare_against.append(QCCompareType.ENRICHMENT_KIT)
        if sequencing_run.enrichment_kit.get_gold_sequencing_runs_qs().exists():
            compare_against.append(QCCompareType.GOLD_ENRICHMENT_KIT_RUNS)

    graph_form = forms.QCCompareTypeForm(columns=SEQUENCING_RUN_QC_COLUMNS,
                                         compare_against=compare_against)
    sample_sheet = sequencing_run.get_current_sample_sheet()
    illumina_qc, show_stats = get_illumina_qc_and_show_stats_for_sample_sheet(sample_sheet)
    read_q30s = None
    try:
        read_q30s = {read: illumina_qc.readq30_set.filter(read=read).get().percent for read in ['R1', 'R2']}
    except:
        pass

    context = {"sequencing_run": sequencing_run,
               'graph_form': graph_form,
               'read_q30s': read_q30s,
               'show_stats': show_stats,
               'illumina_qc': illumina_qc}
    return render(request, 'seqauto/tabs/view_sequencing_run_stats_tab.html', context)


@require_POST
def delete_sequencing_run(request, sequencing_run_id):
    sequencing_run = get_object_or_404(SequencingRun, pk=sequencing_run_id)
    if not request.user.is_superuser:
        raise PermissionDenied()

    sequencing_run.delete()
    return HttpResponse()


@require_POST
def reload_experiment_name(request, sequencing_run_id):
    sequencing_run = get_object_or_404(SequencingRun, pk=sequencing_run_id)

    old_experiment = sequencing_run.experiment

    # TODO: Event log message??
    run_parameters_dir = os.path.join(sequencing_run.path, settings.SEQAUTO_RUN_PARAMETERS_SUB_DIR)
    _, experiment_name = get_run_parameters(run_parameters_dir)
    experiment = None
    if experiment_name:
        experiment, _ = Experiment.objects.get_or_create(name=experiment_name)

    sequencing_run.experiment = experiment
    sequencing_run.save()

    msg = f"{sequencing_run}: experiment {old_experiment} to {experiment}"
    create_event(request.user, "SequencingRun Experiment change", msg)

    return JsonResponse(str(experiment), safe=False)


@require_POST
def assign_data_to_current_sample_sheet(request, sequencing_run_id):
    sequencing_run = get_object_or_404(SequencingRun, pk=sequencing_run_id)
    assign_old_sample_sheet_data_to_current_sample_sheet(request.user, sequencing_run)
    return HttpResponse()


def view_unaligned_reads(request, unaligned_reads_id):
    unaligned_reads = get_object_or_404(UnalignedReads, pk=unaligned_reads_id)

    fastqs = []
    for read_id, fastq in enumerate([unaligned_reads.fastq_r1, unaligned_reads.fastq_r2]):
        try:
            fastqc = fastq.fastqc
        except FastQC.DoesNotExist:
            fastqc = None

        fastqs.append({"read_id": read_id + 1,
                       "fastq": fastq,
                       "fastqc": fastqc})

    context = {"unaligned_reads": unaligned_reads,
               "fastqs": fastqs}
    return render(request, 'seqauto/view_unaligned_reads.html', context)


def view_bam_file(request, bam_file_id):
    bam_file = get_object_or_404(BamFile, pk=bam_file_id)
    form = forms.BamFileForm(instance=bam_file)

    try:
        flagstats = bam_file.flagstats
    except Flagstats.DoesNotExist:
        flagstats = None

    context = {"bam_file": bam_file,
               'form': form,
               'flagstats': flagstats}
    return render(request, 'seqauto/view_bam_file.html', context)


def view_vcf_file(request, vcf_file_id):
    vcf_file = get_object_or_404(VCFFile, pk=vcf_file_id)
    form = forms.VCFFileForm(instance=vcf_file)

    context = {"vcf_file": vcf_file,
               'form': form}
    return render(request, 'seqauto/view_vcf_file.html', context)


def view_combo_vcf_file(request, combo_vcf_file_id):
    combo_vcf_file = get_object_or_404(SampleSheetCombinedVCFFile, pk=combo_vcf_file_id)
    #form = forms.VCFFileForm(instance=vcf_file)

    context = {"combo_vcf_file": combo_vcf_file}
    return render(request, 'seqauto/view_combo_vcf_file.html', context)


def view_qc(request, qc_id):
    qc = get_object_or_404(QC, pk=qc_id)
    form = forms.QCFileForm(instance=qc)

    historical_exec_summaries = list(qc.qcexecsummary_set.all())

    try:
        gene_coverage_collection = qc.qcgenecoverage.gene_coverage_collection
    except QCGeneCoverage.DoesNotExist:
        gene_coverage_collection = None

    context = {"qc": qc,
               'form': form,
               "historical_exec_summaries": historical_exec_summaries,
               "gene_coverage": gene_coverage_collection}
    return render(request, 'seqauto/view_qc.html', context)


def view_qc_exec_summary_tab(request, qc_id):
    qc = get_object_or_404(QC, pk=qc_id)
    graph_form = None
    exec_summary = None
    historical_exec_summaries = list(qc.qcexecsummary_set.all())
    if historical_exec_summaries:
        exec_summary = historical_exec_summaries.pop()
        coverage_columns = list(exec_summary.get_coverage_columns())
        graph_form = forms.QCCompareTypeForm(initial={"compare_against": QCCompareType.SEQUENCING_RUN},
                                             columns=QC_EXEC_SUMMARY_QC_COLUMNS + coverage_columns)

    context = {"qc": qc,
               'graph_form': graph_form,
               "exec_summary": exec_summary}
    return render(request, 'seqauto/tabs/view_qc_exec_summary_tab.html', context)


def view_qc_gene_list_tab(request, qc_id):
    qc = get_object_or_404(QC, pk=qc_id)

    context = {"qc": qc}
    return render(request, 'seqauto/tabs/view_qc_gene_list_tab.html', context)


def view_qc_gene_coverage_collection_tab(request, qc_id):
    qc = get_object_or_404(QC, pk=qc_id)
    try:
        gene_coverage_collection = qc.qcgenecoverage.gene_coverage_collection
    except:
        gene_coverage_collection = None

    context = {"qc": qc,
               "gene_coverage": gene_coverage_collection}
    return render(request, 'seqauto/tabs/view_qc_gene_coverage_collection_tab.html', context)


def view_sample_qc_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    try:
        qc = sample.samplefromsequencingsample.sequencing_sample.get_single_qc()
    except:
        qc = None

    context = {"sample": sample,
               "qc": qc}
    return render(request, 'seqauto/tabs/view_sample_qc_tab.html', context)


def view_enrichment_kit_gene_coverage(request, enrichment_kit_id, gene_symbol):
    user_settings = UserSettings.get_for_user(request.user)
    genome_build = user_settings.default_genome_build

    enrichment_kit = get_object_or_404(EnrichmentKit, pk=enrichment_kit_id)
    context = {"enrichment_kit": enrichment_kit,
               "genome_build": genome_build,
               "gene_symbol": gene_symbol}
    return render(request, 'seqauto/view_enrichment_kit_gene_coverage.html', context)


def view_gold_coverage_summary(request, pk):
    gold_reference = get_object_or_404(GoldReference, pk=pk)
    gold_sequencing_samples_by_run = defaultdict(list)
    qs = gold_reference.goldgenecoveragecollection_set
    for ggcc in qs.order_by(GoldGeneCoverageCollection.SEQUENCING_SAMPLE_PATH):
        gold_sequencing_samples_by_run[ggcc.sequencing_run] = ggcc.sequencing_sample

    context = {"gold_reference": gold_reference,
               "gold_sequencing_samples_by_run": gold_sequencing_samples_by_run}
    return render(request, 'seqauto/view_gold_coverage_summary.html', context)


def enrichment_kits_list(request):
    return render(request, 'seqauto/enrichment_kits_list.html')


def view_enrichment_kit(request, pk):
    enrichment_kit = get_object_or_404(EnrichmentKit, pk=pk)
    default_canonical_transcript_collection = CanonicalTranscriptCollection.get_default()

    context = {"enrichment_kit": enrichment_kit,
               "gold_sequencing_runs_qs": enrichment_kit.get_gold_sequencing_runs_qs(),
               "default_canonical_transcript_collection": default_canonical_transcript_collection}
    return render(request, 'seqauto/view_enrichment_kit.html', context)


def sequencing_stats(request):

    def count_values(klass, column):
        qs = klass.objects.values_list(column).annotate(Count("pk"))
        return list(qs)

    sr_sequencer = count_values(SequencingRun, "sequencer")
    sequencing_run_info = {
        'sequencer': sr_sequencer,
        'sequencer_model': count_values(SequencingRun, "sequencer__sequencer_model"),
    }

    ss_sequencer = count_values(SequencingSample, "sample_sheet__sequencing_run__sequencer")
    sequencing_sample_info = {
        'sequencer': ss_sequencer,
        'sequencer_model':  count_values(SequencingSample, "sample_sheet__sequencing_run__sequencer__sequencer_model"),
        'enrichment_kit': count_values(SequencingSample, "enrichment_kit__name"),
    }

    sample_enrichment_kits_df = get_sample_enrichment_kits_df()
    context = {'num_sequencing_runs': sum([i[1] for i in sr_sequencer]),
               'num_samples': sum([i[1] for i in ss_sequencer]),
               'sequencing_run_info': sequencing_run_info,
               'sequencing_sample_info': sequencing_sample_info,
               'sample_enrichment_kits_df': sample_enrichment_kits_df}
    return render(request, 'seqauto/sequencing_stats.html', context)


def sequencing_stats_data(request):
    context = {"has_sequencing_samples": SequencingSample.objects.all().exists()}
    return render(request, 'seqauto/sequencing_stats_data.html', context)


def qc_data(request):
    return render(request, 'seqauto/qc_data.html')


def sequencing_historical_graphs(request):
    form = forms.QCColumnForm()
    qc_type_totals = dict(QCType.objects.all().values_list("name", "total_field"))

    context = {'form': form,
               'qc_type_totals': qc_type_totals,
               'show_enrichment_kit': QCGraphEnrichmentKitSeparationChoices.SHOW_ENRICHMENT_KIT}
    return render(request, 'seqauto/sequencing_historical_graphs.html', context)


def qc_column_historical_graph(request, qc_column_id, graph_type, enrichment_kit_separation, enrichment_kit_id, use_percent):
    graph_classes = {QCGraphType.LINE_GRAPH: QCColumnLineGraph,
                     QCGraphType.BOX_PLOT: QCColumnBoxPlotGraph}

    graph_class = graph_classes.get(graph_type)
    if not graph_class:
        valid_classes = ','.join(graph_classes.keys())
        msg = f"QCColumn Graph type '{graph_type}' Unknown (should be '{valid_classes}')"
        raise ValueError(msg)
    graph_class_name = full_class_name(graph_class)

    enrichment_kit_separation_dict = dict(QCGraphEnrichmentKitSeparationChoices.choices)
    if enrichment_kit_separation not in enrichment_kit_separation_dict:
        valid_separation = ','.join(enrichment_kit_separation_dict.keys())
        msg = f"QCColumn enrichment_kit_separation '{enrichment_kit_separation}' Unknown (should be '{valid_separation}')"
        raise ValueError(msg)

    if not QCGraphEnrichmentKitSeparationChoices(enrichment_kit_separation).show_enrichment_kit():
        enrichment_kit_id = None

    use_percent = json.loads(use_percent)  # Boolean
    cached_graph = graphcache.async_graph(graph_class_name, qc_column_id, enrichment_kit_separation, enrichment_kit_id, use_percent)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def sequencing_run_qc_graph(request, sequencing_run_id, qc_compare_type):
    _ = QCCompareType(qc_compare_type)  # Check valid
    graph_class_name = full_class_name(SequencingRunQCGraph)
    cached_graph = graphcache.async_graph(graph_class_name, sequencing_run_id, qc_compare_type)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


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

    context = {"current_label": sequencing_run_id,
               "qc_data": sequencing_run_data,
               "label_column": sequencing_run_column,
               "gold_column": gold_column}
    return render(request, 'seqauto/json_graphs/qc_json_graph.html', context)


def index_metrics_qc_graph(request, illumina_qc_id):
    graph_class_name = full_class_name(IndexMetricsQCGraph)
    cached_graph = graphcache.async_graph(graph_class_name, illumina_qc_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def qc_exec_summary_graph(request, qc_exec_summary_id, qc_compare_type):
    graph_class_name = full_class_name(QCExecSummaryGraph)
    cached_graph = graphcache.async_graph(graph_class_name, qc_exec_summary_id, qc_compare_type)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


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
    sample_names = qc_exec_summary_data["qc__bam_file__unaligned_reads__sequencing_sample__sample_name"]
    labels = [get_label(sr, ss) for sr, ss in zip(sequencing_run_names, sample_names)]
    qc_exec_summary_data["label"] = labels

    context = {"current_label": current_label,
               "qc_data": qc_exec_summary_data,
               "label_column": "label",
               "gold_column": gold_column}
    return render(request, 'seqauto/json_graphs/qc_json_graph.html', context)


def get_sequencing_software_versions_template():
    if settings.SEQAUTO_ENABLED:
        base_template = "seqauto/menu_sequencing_data_base.html"
    else:
        base_template = "snpdb/menu/menu_settings_base.html"
    return base_template


def sequencing_software_versions(request):
    # TODO: Forms etc
    context = {"base_template": get_sequencing_software_versions_template()}
    return render(request, 'seqauto/sequencing_software_versions.html', context)


@method_decorator([staff_only], name='dispatch')
class SeqautoUpdateView(UpdateView):
    template_name = "snpdb/update_form.html"
    widgets = {}

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = get_sequencing_software_versions_template()
        context['title'] = str(self.model.__name__)
        return context

    def get_form(self, form_class=None):
        form = super().get_form(form_class)
        for f in form.fields:
            w = self.widgets.get(f)
            if w:
                form.fields[f].widget = w
        return form


class LibraryUpdate(SeqautoUpdateView):
    model = Library
    fields = get_model_fields(Library)
    #widgets = {"name" : TextInput()}


class SequencerUpdate(SeqautoUpdateView):
    model = Sequencer
    fields = get_model_fields(Sequencer)


class AssayUpdate(SeqautoUpdateView):
    model = Assay
    fields = get_model_fields(Assay)


class AlignerUpdate(SeqautoUpdateView):
    model = Aligner
    fields = get_model_fields(Aligner)


class VariantCallerUpdate(SeqautoUpdateView):
    model = VariantCaller
    fields = get_model_fields(VariantCaller)


class VariantCallingPipelineUpdate(SeqautoUpdateView):
    model = VariantCallingPipeline
    fields = get_model_fields(VariantCallingPipeline)
