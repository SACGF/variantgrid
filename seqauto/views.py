import os
from collections import defaultdict

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.db.models.query_utils import Q
from django.http.response import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from eventlog.models import create_event
from genes.models import CanonicalTranscriptCollection
from library.log_utils import log_traceback
from seqauto import forms
from seqauto.forms import AllEnrichmentKitForm, AutocompleteSequencingRunForm, SequencingRunForm
from seqauto.illumina.run_parameters import get_run_parameters
from seqauto.models import (
    QC,
    BamFile,
    EnrichmentKit,
    Experiment,
    FastQC,
    Flagstats,
    GoldGeneCoverageCollection,
    GoldReference,
    JointCalledVCF,
    QCGeneCoverage,
    QCType,
    SequencingRun,
    SingleSampleVCF,
    UnalignedReads,
)
from seqauto.models.models_enums import QCCompareType
from seqauto.qc.sequencing_run_utils import SEQUENCING_RUN_QC_COLUMNS
from seqauto.sequencing_files.sample_sheet import (
    assign_old_sample_sheet_data_to_current_sample_sheet,
)
from snpdb.models import Sample, UserSettings


def sequencing_data(request):
    return render(request, 'seqauto/sequencing_data.html')


def experiments(request):
    return render(request, 'seqauto/experiments.html')


def sequencing_runs(request):
    context = {
        "enrichment_kit_form": AllEnrichmentKitForm(),
        "sequencing_run_form": AutocompleteSequencingRunForm(),
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


def get_illumina_qc_and_show_stats_for_sample_sheet(sample_sheet):
    try:
        illumina_qc = sample_sheet.illuminaflowcellqc
        show_stats = True
    except Exception:
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

    vcf_types = {
        "Joint Called VCF": Q(vcf__uploadedvcf__backendvcf__joint_called_vcf__isnull=False),
        "Single Sample VCF": Q(vcf__uploadedvcf__backendvcf__single_sample_vcf__isnull=False),
    }

    run_vcfs = defaultdict(list)
    sr_vcf_qs = sequencing_run.vcffromsequencingrun_set.all()
    for vcf_type, q in vcf_types.items():
        for vcf_for_run in sr_vcf_qs.filter(q):
            run_vcfs[vcf_type].append((vcf_for_run.get_variant_caller(), vcf_for_run.vcf, vcf_for_run.vcf.can_view(request.user)))

    context = {
        "sequencing_run": sequencing_run,
        "sequencing_run_form": sequencing_run_form,
        'tab_id': tab_id,
        'run_vcfs': dict(run_vcfs),
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
    except Exception:
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
        read_q30s = {read: illumina_qc.readq30_set.get(read=read).percent for read in ['R1', 'R2']}
    except Exception:
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


def view_single_sample_vcf(request, single_sample_vcf_id):
    single_sample_vcf = get_object_or_404(SingleSampleVCF, pk=single_sample_vcf_id)
    context = {"single_sample_vcf": single_sample_vcf}
    return render(request, 'seqauto/view_single_sample_vcf.html', context)


def view_vcf_file(request, vcf_file_id):
    return view_single_sample_vcf(request, vcf_file_id)


def view_joint_called_vcf(request, joint_called_vcf_id):
    joint_called_vcf = get_object_or_404(JointCalledVCF, pk=joint_called_vcf_id)
    context = {"joint_called_vcf": joint_called_vcf}
    return render(request, 'seqauto/view_joint_called_vcf.html', context)


def view_combo_vcf_file(request, combo_vcf_file_id):
    return view_joint_called_vcf(request, combo_vcf_file_id)


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
        exec_summary_qc = QCType.objects.get(name="ExecSummaryQC")
        qc_exec_summary_columns = list(exec_summary_qc.qccolumn_set.all().values_list("field", flat=True))

        exec_summary = historical_exec_summaries.pop()
        coverage_columns = exec_summary.get_coverage_columns()
        graph_form = forms.QCCompareTypeForm(initial={"compare_against": QCCompareType.SEQUENCING_RUN},
                                             columns=qc_exec_summary_columns + coverage_columns)

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
    except Exception:
        gene_coverage_collection = None

    context = {"qc": qc,
               "gene_coverage": gene_coverage_collection}
    return render(request, 'seqauto/tabs/view_qc_gene_coverage_collection_tab.html', context)


def view_sample_qc_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    try:
        qc = sample.samplefromsequencingsample.sequencing_sample.get_single_qc()
    except Exception:
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
