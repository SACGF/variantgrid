import logging
from collections import Counter, defaultdict

from django.contrib import messages
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_POST

from annotation.models import (
    AnnotationRun,
    InvalidAnnotationVersionError,
    VariantAnnotationVersion,
)
from annotation.models.models import AnnotationPipelineVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.pipelines import get_pipeline, get_runner, versioned_pipeline_types
from annotation.tasks.annotate_variants import annotation_run_retry
from annotation.tasks.annotation_scheduler_task import (
    annotation_scheduler,
    dispatch_annotation_runs,
    subdivide_annotation_range_lock,
)
from library.django_utils import get_field_counts, require_superuser
from snpdb.models import GenomeBuild, JobsControl
from variantgrid.celery import app
from variantgrid.deployment_validation.annotation_status_checks import (
    get_variant_annotation_progress,
    is_variant_annotation_version_populated,
)


def _register_pipeline_version_post_name(pipeline_type, genome_build) -> str:
    return f"register-pipeline-version-{pipeline_type}-{genome_build.name}"


def _handle_register_pipeline_version_post(request):
    """ #720: the runs page equivalent of create_new_annotation_pipeline_version - reads the installed
        tool's version and registers it, so rolling out an upgrade doesn't need shell access. """
    for pipeline_type in versioned_pipeline_types():
        runner = get_runner(pipeline_type)
        label = VariantAnnotationPipelineType(pipeline_type).label
        for genome_build in GenomeBuild.builds_with_annotation():
            if _register_pipeline_version_post_name(pipeline_type, genome_build) not in request.POST:
                continue
            if not get_pipeline(pipeline_type).enabled:
                messages.add_message(request, messages.ERROR, f"{label} is not enabled in settings")
                continue
            if not runner.supports_genome_build(genome_build):
                messages.add_message(request, messages.ERROR,
                                     f"{label} does not support {genome_build}")
                continue
            try:
                pipeline_version, created = runner.get_or_create_current_version(genome_build)
            except OSError as e:
                # Reading the version shells out to the tool, which the web worker may not have installed
                messages.add_message(request, messages.ERROR,
                                     f"{label} ({genome_build}) - couldn't read the installed version: {e}")
                continue

            if not created:
                messages.add_message(request, messages.INFO,
                                     f"{label} ({genome_build}) - installed tool is already registered as "
                                     f"{pipeline_version} [{pipeline_version.get_status_display()}]")
            elif pipeline_version.is_active:
                annotation_scheduler.si().apply_async()
                messages.add_message(request, messages.INFO,
                                     f"{label} ({genome_build}) - registered {pipeline_version} as ACTIVE "
                                     f"and queued the scheduler to create its annotation runs")
            else:
                messages.add_message(request, messages.INFO,
                                     f"{label} ({genome_build}) - registered {pipeline_version} as NEW. "
                                     f"Promote it to re-annotate every range with it.")


@require_superuser
def variant_annotation_runs(request):
    as_display = dict(AnnotationStatus.choices)

    genome_build_field_counts = defaultdict(dict)
    genome_build_summary = defaultdict(dict)
    genome_build_summary_combined = defaultdict(dict)

    if request.method == "POST":
        if "unpause-jobs" in request.POST:
            JobsControl.resume(by=str(request.user))
            messages.add_message(request, messages.INFO,
                                 "Resumed analysis + annotation job dispatch")
        if "pause-jobs" in request.POST:
            JobsControl.pause(reason=f"Paused from annotation runs page by {request.user}",
                              by=str(request.user))
            messages.add_message(request, messages.WARNING,
                                 "Paused analysis + annotation job dispatch")

        if "annotation-scheduler" in request.POST:
            annotation_scheduler.si().apply_async()

        if "run-scheduler-new" in request.POST:
            annotation_scheduler.si(status=VariantAnnotationVersion.Status.NEW).apply_async()
            messages.add_message(request, messages.INFO, "Queued annotation scheduler against NEW")

        # #720: promoting a NEW AnnotationPipelineVersion is what starts a supplementary pipeline's
        # re-annotation - the scheduler then creates a run per range lock against it, the same query that
        # backfills a newly-enabled pipeline.
        for pipeline_version in AnnotationPipelineVersion.objects.filter(
                status=AnnotationPipelineVersion.Status.NEW):
            if f"promote-pipeline-version-{pipeline_version.pk}" in request.POST:
                pipeline_version.promote_to_active()
                annotation_scheduler.si().apply_async()
                messages.add_message(request, messages.INFO,
                                     f"Promoted {pipeline_version} to ACTIVE - queued the scheduler to "
                                     f"create its annotation runs")

        _handle_register_pipeline_version_post(request)

        for genome_build in GenomeBuild.builds_with_annotation():
            if f"promote-to-active-{genome_build.name}" in request.POST:
                new_vav = VariantAnnotationVersion.latest(genome_build,
                                                         status=VariantAnnotationVersion.Status.NEW)
                if new_vav is None:
                    messages.add_message(request, messages.ERROR,
                                         f"{genome_build} - no NEW VariantAnnotationVersion to promote")
                elif not is_variant_annotation_version_populated(new_vav):
                    messages.add_message(request, messages.ERROR,
                                         f"{genome_build} - NEW VAV pk={new_vav.pk} is not yet populated")
                else:
                    try:
                        new_vav.promote_to_active()
                        messages.add_message(request, messages.INFO,
                                             f"{genome_build} - promoted VAV pk={new_vav.pk} to ACTIVE")
                    except (ValueError, InvalidAnnotationVersionError) as e:
                        messages.add_message(request, messages.ERROR, str(e))

            for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build):
                annotation_runs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
                message = None
                if f"set-non-finished-to-error-{genome_build.name}-{vav.pk}" in request.POST:
                    num_errored = 0
                    non_finished_statuses = [AnnotationStatus.FINISHED, AnnotationStatus.ERROR]
                    for annotation_run in annotation_runs.exclude(status__in=non_finished_statuses):
                        if celery_task := annotation_run.task_id:
                            logging.info("Terminating celery job '%s'", celery_task)
                            app.control.revoke(celery_task, terminate=True)  # @UndefinedVariable
                        annotation_run.error_exception = "Manually failed"
                        annotation_run.save()
                        num_errored += 1
                    message = f"{genome_build} - set {num_errored} annotation runs to Error"
                elif f"retry-annotation-runs-{genome_build.name}-{vav.pk}" in request.POST:
                    num_retrying = 0
                    for annotation_run in annotation_runs.filter(status=AnnotationStatus.ERROR):
                        annotation_run_retry(annotation_run)
                        num_retrying += 1
                    message = f"{genome_build} - retrying {num_retrying} annotation runs."

                if message:
                    messages.add_message(request, messages.INFO, message)

    for genome_build in GenomeBuild.builds_with_annotation():
        active_pipeline_versions = {
            pipeline_type: AnnotationPipelineVersion.get_active(pipeline_type, genome_build)
            for pipeline_type in versioned_pipeline_types()
        }
        for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build).order_by("-annotation_date"):
            base_qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
            # Split by pipeline type so standard-variant work (the bottleneck) is shown separately from
            # structural-variant runs (mostly empty/instant). Only include types that have runs.
            summary_by_type = {}
            field_counts_by_type = {}
            combined_summary = Counter()  # across pipeline types - drives the per-VAV action buttons
            for pipeline_type in VariantAnnotationPipelineType:
                type_qs = base_qs.filter(pipeline_type=pipeline_type)
                # #720: a versioned pipeline keeps the runs of every tool version it has been through, so
                # scope the outstanding/finished counts to the one being scheduled against now. The
                # earlier versions' runs are history, not work.
                if active_version := active_pipeline_versions.get(pipeline_type):
                    type_qs = type_qs.filter(pipeline_version=active_version)
                field_counts = get_field_counts(type_qs, "status")
                if not field_counts:
                    continue
                summary_data = Counter()
                for field, count in field_counts.items():
                    summary = AnnotationStatus.get_summary_state(field)
                    summary_data[summary] += count
                    combined_summary[summary] += count
                summary_by_type[pipeline_type.value] = summary_data
                field_counts_by_type[pipeline_type.value] = {as_display[k]: v for k, v in field_counts.items()}

            genome_build_summary[genome_build.pk][vav.pk] = summary_by_type
            genome_build_summary_combined[genome_build.pk][vav.pk] = combined_summary
            genome_build_field_counts[genome_build.pk][vav] = field_counts_by_type

    current_variant_annotation_versions = VariantAnnotationVersion.latest_for_all_builds()
    new_variant_annotation_versions = VariantAnnotationVersion.objects.filter(
        status=VariantAnnotationVersion.Status.NEW).order_by("genome_build__name", "annotation_date")
    historical_variant_annotation_versions = VariantAnnotationVersion.objects.filter(
        status=VariantAnnotationVersion.Status.HISTORICAL).order_by("genome_build__name", "-annotation_date")

    genome_build_status_panel = {}
    for genome_build in GenomeBuild.builds_with_annotation():
        active_vav = VariantAnnotationVersion.latest(genome_build,
                                                    status=VariantAnnotationVersion.Status.ACTIVE)
        new_vav = VariantAnnotationVersion.latest(genome_build,
                                                 status=VariantAnnotationVersion.Status.NEW)
        new_progress = None
        new_populated = False
        new_gene_annotation_blocker = None
        if new_vav is not None:
            new_progress = get_variant_annotation_progress(new_vav)
            new_populated = is_variant_annotation_version_populated(new_vav)
            new_gene_annotation_blocker = new_vav.get_gene_annotation_promote_blocker()
        historical_count = VariantAnnotationVersion.objects.filter(
            genome_build=genome_build, status=VariantAnnotationVersion.Status.HISTORICAL
        ).count()
        genome_build_status_panel[genome_build.name] = {
            "active": active_vav,
            "new": new_vav,
            "new_progress": new_progress,
            "new_populated": new_populated,
            "new_gene_annotation_blocker": new_gene_annotation_blocker,
            "historical_count": historical_count,
        }

    # #720: the AnnotationPipelineVersion equivalent of the VAV status panel - what each supplementary
    # pipeline is annotating against, and whether an upgrade is registered and waiting to be promoted.
    pipeline_version_panels = []
    for genome_build in GenomeBuild.builds_with_annotation():
        for pipeline_type in versioned_pipeline_types():
            pipeline_version_panels.append({
                "genome_build": genome_build,
                "label": VariantAnnotationPipelineType(pipeline_type).label,
                "enabled": get_pipeline(pipeline_type).enabled,
                "supported": get_runner(pipeline_type).supports_genome_build(genome_build),
                "active": AnnotationPipelineVersion.get_active(pipeline_type, genome_build),
                "new": AnnotationPipelineVersion.get_new(pipeline_type, genome_build),
                "historical_count": AnnotationPipelineVersion.objects.filter(
                    pipeline_type=pipeline_type, genome_build=genome_build,
                    status=AnnotationPipelineVersion.Status.HISTORICAL).count(),
                "register_name": _register_pipeline_version_post_name(pipeline_type, genome_build),
                "register_command": f"python3 manage.py create_new_annotation_pipeline_version "
                                    f"--pipeline-type={pipeline_type} --genome-build={genome_build}",
            })

    any_new_vav = VariantAnnotationVersion.objects.filter(
        status=VariantAnnotationVersion.Status.NEW
    ).exists()

    # Don't force-create the singleton just by viewing the page
    jobs_control = JobsControl.objects.filter(pk=JobsControl.SINGLETON_PK).first()

    context = {
        "genome_build_summary": dict(genome_build_summary),
        "genome_build_summary_combined": dict(genome_build_summary_combined),
        "genome_build_field_counts": dict(genome_build_field_counts),
        "pipeline_type_labels": dict(VariantAnnotationPipelineType.choices),
        "current_variant_annotation_versions": current_variant_annotation_versions,
        "new_variant_annotation_versions": new_variant_annotation_versions,
        "historical_variant_annotation_versions": historical_variant_annotation_versions,
        "genome_build_status_panel": genome_build_status_panel,
        "pipeline_version_panels": pipeline_version_panels,
        "any_new_vav": any_new_vav,
        "jobs_control": jobs_control,
        "jobs_paused": bool(jobs_control and jobs_control.paused),
    }
    return render(request, "annotation/variant_annotation_runs.html", context)


@require_superuser
def view_annotation_run(request, annotation_run_id):
    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)

    # There may be other runs of different types (don't care about them)
    other_annotation_runs_qs = AnnotationRun.objects.filter(
        annotation_range_lock=annotation_run.annotation_range_lock,
        pipeline_type=annotation_run.pipeline_type).exclude(pk=annotation_run.pk)
    can_retry_annotation_run = not other_annotation_runs_qs.exists()
    can_retry_annotation_run_upload = can_retry_annotation_run and annotation_run.vcf_annotated_filename
    can_subdivide_annotation_run = False
    if arl := annotation_run.annotation_range_lock:
        can_subdivide_annotation_run = arl.can_subdivide()

    context = {
        "annotation_run": annotation_run,
        "can_retry_annotation_run": can_retry_annotation_run,
        "can_retry_annotation_run_upload": can_retry_annotation_run_upload,
        "can_subdivide_annotation_run": can_subdivide_annotation_run,
        "can_make_annotation_run_local": annotation_run.get_status() == AnnotationStatus.EXTERNAL_DUMP_COMPLETED,
    }
    return render(request, "annotation/view_annotation_run.html", context)


@require_superuser
@require_POST
def retry_annotation_run(request, annotation_run_id, upload_only=False):
    """ Deletes - then re-tries using annotation lock """
    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)
    annotation_run = annotation_run_retry(annotation_run, upload_only=upload_only)

    msg = 'Re-trying annotation'
    if upload_only:
        msg += " (upload only)"
    status = messages.INFO
    messages.add_message(request, status, msg, extra_tags='import-message')
    return redirect(annotation_run)


@require_superuser
@require_POST
def retry_annotation_run_upload(request, annotation_run_id):
    return retry_annotation_run(request, annotation_run_id, upload_only=True)


@require_superuser
@require_POST
def make_annotation_run_local(request, annotation_run_id):
    """ #1568: revert an external annotation run back to the local pipeline so it is annotated
        locally, regardless of size. Clears external + dump state (-> CREATED) and kicks the dispatcher. """
    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)
    if annotation_run.get_status() == AnnotationStatus.EXTERNAL_DUMP_COMPLETED:
        annotation_run.revert_external_to_local()
        dispatch_annotation_runs.si(annotation_run.variant_annotation_version.pk).apply_async()
        messages.add_message(request, messages.INFO,
                             "Reverted to local pipeline - will be annotated locally",
                             extra_tags='import-message')
    else:
        messages.add_message(request, messages.WARNING,
                             "Annotation run is not awaiting external annotation - nothing to do",
                             extra_tags='import-message')
    return redirect(annotation_run)


@require_superuser
@require_POST
def subdivide_annotation_run(request, annotation_run_id):
    """ Sometimes runs fail w/out of memory etc (perhaps due to too many transcripts) - be able to subdivide """

    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)
    pipeline_type = annotation_run.pipeline_type
    pipeline_version = annotation_run.pipeline_version  # #720: the halves stay on the version this run had
    old_range_lock = annotation_run.annotation_range_lock
    version_id = old_range_lock.version_id

    # subdivide_annotation_range_lock deletes the run(s) on this lock, shrinks it to the bottom half and
    # returns a new lock for the top half. Create a fresh CREATED run per half (each with its lock set,
    # so no rangeless run is ever committed - #1654) and hand back to the dispatcher.
    new_range_lock = subdivide_annotation_range_lock(old_range_lock)
    bottom_half_run = AnnotationRun.objects.create(annotation_range_lock=old_range_lock,
                                                   pipeline_type=pipeline_type,
                                                   pipeline_version=pipeline_version)
    AnnotationRun.objects.create(annotation_range_lock=new_range_lock, pipeline_type=pipeline_type,
                                 pipeline_version=pipeline_version)
    dispatch_annotation_runs.si(version_id).apply_async()

    # The original run was deleted by the subdivision - land on the bottom-half replacement.
    return redirect(bottom_half_run)
