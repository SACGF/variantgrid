
from django.contrib import messages
from django.shortcuts import render

from library.django_utils import (
    get_field_counts,
    require_superuser,
)
from snpdb.models import (
    Allele,
    AlleleLiftover,
    GenomeBuild,
    LiftoverRun,
)
from snpdb.models.models_enums import (
    AlleleConversionTool,
    ProcessingStatus,
)
from snpdb.tasks.liftover_tasks import liftover_alleles


@require_superuser
def liftover_runs(request):
    genome_builds = GenomeBuild.builds_with_annotation()
    if request.method == 'POST':
        dest_genome_build = None
        for genome_build in genome_builds:
            if f"liftover_to_{genome_build.name}" in request.POST:
                dest_genome_build = genome_build
                break
        if dest_genome_build is None:
            raise ValueError("Could not determine dest genome build from liftover_runs POST")
        messages.add_message(request, messages.INFO, f"Lifting over alleles to {dest_genome_build}")
        liftover_alleles.si(request.user.username, dest_genome_build.name).apply_async()

    alleles_missing_variants = {}
    for genome_build in genome_builds:
        alleles_missing_variants[genome_build.name] = Allele.missing_variants_for_build(genome_build).count()

    qs_allele_liftover = AlleleLiftover.objects.all()

    tool_status = {}
    failure_rate = {}
    for act in AlleleConversionTool:
        qs_act = qs_allele_liftover.filter(liftover__conversion_tool=act)
        if act_data := get_field_counts(qs_act, "status"):
            tool_status[act.label] = {ProcessingStatus(k): v for k, v in act_data.items()}
            num_error = act_data.get(ProcessingStatus.ERROR, 0)
            num_success = act_data.get(ProcessingStatus.SUCCESS, 0)
            total = num_error + num_success
            if total:
                failure_rate[act.label] = 100.0 * num_error / total

    used_states = set()
    for s in [d.keys() for d in tool_status.values()]:
        used_states.update(s)

    processing_status_cols = [ps for ps in ProcessingStatus if ps in used_states]  # In ProcessingStatus order
    context = {
        "genome_builds": genome_builds,
        "alleles_missing_variants": alleles_missing_variants,
        "processing_status_cols": processing_status_cols,
        "tool_status": tool_status,
        "failure_rate": failure_rate,
    }

    return render(request, "snpdb/liftover/liftover_runs.html", context)


@require_superuser
def view_liftover_run(request, liftover_run_id):
    liftover_run = LiftoverRun.objects.get(pk=liftover_run_id)
    qs = liftover_run.alleleliftover_set.all()
    raw_status_counts = get_field_counts(qs, "status")
    status_counts = {ProcessingStatus(k).label: v for k, v in raw_status_counts.items() if k is not None}

    context = {
        "liftover_run": liftover_run,
        "status_counts": status_counts,
    }
    return render(request, "snpdb/liftover/view_liftover_run.html", context)
