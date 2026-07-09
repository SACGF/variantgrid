import gzip
import logging
import os

import celery
from celery import chain
from celery.canvas import Signature
from django.conf import settings
from django.db.models.functions.math import Abs
from django.db.models.query_utils import Q
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.annotsv_annotation import run_annotsv, annotsv_check_command_line_version_match
from annotation.models import AnnotationStatus, GenomeBuild, VariantAnnotationPipelineType
from annotation.models.models import AnnotationRun, InvalidAnnotationVersionError
from annotation.signals.manual_signals import annotation_run_complete_signal
from annotation.vcf_files.import_vcf_annotations import import_vcf_annotations
from annotation.vep_annotation import get_vep_command, vep_check_command_line_version_match
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.log_utils import get_traceback, report_message, log_traceback
from library.utils import execute_cmd
from library.utils.file_utils import name_from_filename, mk_path_for_file
from snpdb.variants_to_vcf import write_contig_sorted_values_to_vcf_file, VARIANT_GRID_INFO_DICT


# #2667: kick the single-authority dispatcher by name to avoid importing annotation_scheduler_task
# (which imports this module). Mirror of analysis _trigger_rescheduling (#346).
DISPATCH_ANNOTATION_RUNS_TASK = "annotation.tasks.annotation_scheduler_task.dispatch_annotation_runs"


def _trigger_dispatch(variant_annotation_version_id):
    """ A run just completed (success OR failure) - a worker has freed up. Kick the dispatcher to
        launch the next (merged) batch. Both kicks serialise through scheduling_single_worker and
        fast-exit if nothing is dispatchable, so firing twice is safe (the dispatcher is single
        authority and leases atomically).

        Why two kicks: race condition releasing the lock. This worker has just cleared the run's
        task_id/lease and saved, but the dispatcher (a separate worker) computes free capacity by
        reading in-flight runs from the DB. If its read transaction starts before our completion
        commit is visible to it, it still counts this run as in-flight, sees no free slot, and exits -
        leaving the freed slot idle until the next event. The immediate kick handles the common case
        (commit already visible, lowest latency); the 3s-delayed kick re-runs the dispatcher once the
        lock release is unambiguously committed, so the freed capacity is actually picked up. """
    sig = Signature(DISPATCH_ANNOTATION_RUNS_TASK, args=(variant_annotation_version_id,))
    sig.apply_async()
    sig.apply_async(countdown=3)


def _dispatch_trigger_sig(variant_annotation_version_id) -> Signature:
    """ Immutable by-name dispatcher kick, for chaining after a retry's cleanup task so the dispatcher
        only re-picks the run once the cleanup has committed (a retry must not race its own scrub). The
        run then routes by status - upload-only stays ANNOTATION_COMPLETED (import lane), a full retry is
        a fresh CREATED run (VEP lane). #1649 """
    return Signature(DISPATCH_ANNOTATION_RUNS_TASK, args=(variant_annotation_version_id,), immutable=True)


@celery.shared_task
def delete_annotation_run(annotation_run_id):
    try:
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        annotation_run.status = AnnotationStatus.DELETING
        annotation_run.save()
        annotation_run.delete()
    except:
        log_traceback()
        raise


@celery.shared_task
def delete_annotation_run_uploaded_data(annotation_run_id):
    """ Deletes related objects but not actual run (used for retry annotation upload) """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    annotation_run.delete_related_objects()


@celery.shared_task
def assign_range_lock_to_annotation_run(annotation_run_id, annotation_range_lock_id):
    # has 1-to-1 so any previous AnnotationRuns linking to AnnotationRangeLock must have been deleted
    AnnotationRun.objects.filter(pk=annotation_run_id).update(annotation_range_lock_id=annotation_range_lock_id)


@celery.shared_task
def annotate_variants(annotation_run_id):
    """ VEP lane (annotation_workers): dump unannotated variants + run the VEP subprocess, then stop at
        ANNOTATION_COMPLETED and release the lease. The DB upload is deliberately NOT done here - the
        dispatcher re-picks the completed run in its resume lane and runs import_annotation_run on
        db_workers, so a throttled VEP slot is never held through the bulk insert. External runs
        (#1568) never reach this - they're dumped off-VM and rejoin post-VEP as ANNOTATION_COMPLETED.
        See #1649. """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
    logging.info("annotate_variants: %s", annotation_run)

    # External annotation (#1568): VEP for these runs is managed off-VM via the annotation_external
    # command. Never auto-run VEP here while waiting for the operator to import an annotated VCF.
    if annotation_run.external and annotation_run.vcf_annotated_filename is None:
        logging.info("Skipping external AnnotationRun %s (awaiting external annotation)", annotation_run.pk)
        return

    # task_id used as Celery lock
    num_modified = AnnotationRun.objects.filter(pk=annotation_run.pk,
                                                task_id__isnull=True).update(task_id=annotate_variants.request.id)
    if num_modified != 1:
        msg = f"Celery couldn't get task_id lock on AnnotationRun: {annotation_run.pk}"
        annotation_run.celery_task_logs[annotate_variants.request.id] = msg
        annotation_run.save()
        raise ValueError(msg)

    try:
        # Reload to get updated task_id
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        if blocker := annotation_run.variant_annotation_version.get_annotation_run_blocker():
            # We need this so that transcript/versions are in DB so FKs link
            msg = f"{annotation_run.variant_annotation_version} {blocker}"
            raise InvalidAnnotationVersionError(msg)
        annotation_run.task_id = annotate_variants.request.id
        annotation_run.set_task_log("start", timezone.now())
        annotation_run.save()

        if annotation_run.vcf_annotated_filename is None:
            dump_and_annotate_variants(annotation_run)
        # DB upload now runs as a separate db_workers task (import_annotation_run), launched by the
        # dispatcher when this run reaches ANNOTATION_COMPLETED. #1649
    except Exception as e:
        tb = get_traceback()
        error_message = f"{e}: {annotation_run.pipeline_stderr}"

        name = 'Annotation pipeline run ' + str(annotation_run.id)
        report_message(message=name,
                       level='error',
                       extra_data={'output': error_message,
                                   'error': tb})

        create_event(None, "AnnotationRun failed", tb, severity=LogLevel.ERROR)
        annotation_run.error_exception = tb
        annotation_run.set_task_log("error_exception", tb)
        annotation_run.save()
        raise
    finally:
        annotation_run.set_task_log("end", timezone.now())
        # #2667: release the dispatcher lease so the now-completed run no longer looks in-flight,
        # then kick the dispatcher to refill this freed worker slot (covers success and failure).
        annotation_run.task_id = None
        annotation_run.leased_by = None
        annotation_run.lease_expires = None
        annotation_run.save()
        _trigger_dispatch(annotation_run.annotation_range_lock.version_id)


@celery.shared_task
def import_annotation_run(annotation_run_id):
    """ Import lane (db_workers): bulk-load an already-annotated VCF into the DB. Reached for every run
        once it is ANNOTATION_COMPLETED with an annotated VCF present - whether VEP just finished in-VM,
        an external run (#1568) was imported, or an upload-only retry reset it. Kept off annotation_workers
        so quick DB inserts never consume a throttled VEP slot. See #1649. """
    annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)

    # task_id used as Celery lock (as annotate_variants) so a double-dispatch can't double-import
    num_modified = AnnotationRun.objects.filter(
        pk=annotation_run.pk, task_id__isnull=True).update(task_id=import_annotation_run.request.id)
    if num_modified != 1:
        msg = f"Celery couldn't get task_id lock on AnnotationRun: {annotation_run.pk}"
        annotation_run.celery_task_logs[import_annotation_run.request.id] = msg
        annotation_run.save()
        raise ValueError(msg)

    try:
        # Reload to get updated task_id
        annotation_run = AnnotationRun.objects.get(pk=annotation_run_id)
        annotation_run.task_id = import_annotation_run.request.id
        annotation_run.set_task_log("import_start", timezone.now())
        annotation_run.save()

        import_vcf_annotations(annotation_run)
        # The run is only truly complete once imported (moved here from annotate_variants). #1649
        annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                            variant_annotation_version=annotation_run.annotation_range_lock.version)
    except Exception:
        tb = get_traceback()
        create_event(None, "AnnotationRun import failed", tb, severity=LogLevel.ERROR)
        annotation_run.error_exception = tb
        annotation_run.set_task_log("error_exception", tb)
        annotation_run.save()
        raise
    finally:
        annotation_run.set_task_log("import_end", timezone.now())
        # #2667: release the dispatcher lease so the now-completed run no longer looks in-flight, then
        # kick the dispatcher to refill this freed import slot (covers success and failure).
        annotation_run.task_id = None
        annotation_run.leased_by = None
        annotation_run.lease_expires = None
        annotation_run.save()
        _trigger_dispatch(annotation_run.annotation_range_lock.version_id)


def dump_variants(annotation_run, dump_dir=None) -> int:
    """ Write the unannotated variants in range to the dump VCF and set dump_* fields; returns dump count.
        Factored out of dump_and_annotate_variants so the annotation_external --dump command (#1568) can
        dump (into --output-dir) and stop before VEP. """
    vcf_dump_filename = annotation_run.get_dump_filename(dump_dir=dump_dir)
    annotation_run.dump_start = timezone.now()
    annotation_run.vcf_dump_filename = vcf_dump_filename
    annotation_run.save()

    genome_build = annotation_run.genome_build
    vcf_dump_count = _unannotated_variants_to_vcf(genome_build, vcf_dump_filename,
                                                  annotation_run.annotation_range_lock,
                                                  annotation_run.pipeline_type)

    annotation_run.dump_count = vcf_dump_count
    annotation_run.dump_end = timezone.now()
    annotation_run.save()
    return vcf_dump_count


def dump_and_annotate_variants(annotation_run, vep_version_check=True):
    if vep_version_check:
        # Do a check before we annotate
        vep_check_command_line_version_match(annotation_run.variant_annotation_version)
        # Counterpart for AnnotSV - skipped entirely when not enabled
        if (settings.ANNOTATION_ANNOTSV_ENABLED
                and annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
            annotsv_check_command_line_version_match(annotation_run.variant_annotation_version)

    vcf_dump_count = dump_variants(annotation_run)
    vcf_dump_filename = annotation_run.vcf_dump_filename

    genome_build = annotation_run.genome_build
    annotation_consortium = annotation_run.annotation_consortium

    logging.info("Annotating %d variants", vcf_dump_count)
    if vcf_dump_count:
        name = name_from_filename(vcf_dump_filename)
        vcf_annotated_basename = f"{name}.vep_annotated_{genome_build.name}.vcf.gz"
        vcf_annotated_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, vcf_annotated_basename)

        cmd = get_vep_command(vcf_dump_filename, vcf_annotated_filename, genome_build, annotation_consortium,
                              annotation_run.pipeline_type,
                              variant_annotation_version=annotation_run.variant_annotation_version)
        annotation_run.annotation_start = timezone.now()
        annotation_run.pipeline_command = " ".join(cmd)
        annotation_run.save()
        return_code, std_out, std_err = execute_cmd(cmd)
        # VEP can produce enormous output (>1GB) for large batches - PostgreSQL has a 1GB field limit
        max_output = 1_000_000
        annotation_run.pipeline_stdout = std_out[:max_output] if std_out else std_out
        annotation_run.pipeline_stderr = std_err[:max_output] if std_err else std_err
        logging.info(f"VEP returned code: {return_code}")

        if return_code != 0:
            annotation_run.save()  # save stdout/stderr
            raise RuntimeError(f"VEP returned {return_code}")

        annotation_run.vcf_annotated_filename = vcf_annotated_filename
        annotation_run.annotation_end = timezone.now()

        if (settings.ANNOTATION_ANNOTSV_ENABLED
                and annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
            annotsv_dir = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR,
                                       f"annotsv_{annotation_run.pk}")
            tsv, rc, stdout, stderr = run_annotsv(vcf_dump_filename, annotsv_dir,
                                                  genome_build, annotation_consortium)
            if rc == 0 and os.path.exists(tsv):
                annotation_run.annotsv_tsv_filename = tsv
                annotation_run.annotsv_error = None
            else:
                tsv_missing = "" if os.path.exists(tsv) else f" (expected TSV not found: {tsv})"
                error_blob = (
                    f"rc={rc}{tsv_missing}\n"
                    f"--- stderr ---\n{stderr or ''}\n"
                    f"--- stdout ---\n{stdout or ''}"
                )
                annotation_run.annotsv_error = error_blob[:100_000]
                logging.warning(
                    "AnnotSV stage failed for AnnotationRun %s: rc=%s%s\nstderr:\n%s\nstdout:\n%s",
                    annotation_run.pk, rc, tsv_missing,
                    (stderr or "")[-4000:], (stdout or "")[-4000:],
                )
    else:
        # Now we have standard/CNV type pipelines, it's possible some can be empty
        annotation_run.annotated_count = 0
        annotation_run.annotation_end = timezone.now()

    annotation_run.save()


def annotation_run_retry(annotation_run: AnnotationRun, upload_only=False) -> AnnotationRun:
    if upload_only and annotation_run.vcf_annotated_filename is None:
        msg = "Retry annotation run upload only requires annotation VCF to be written"
        raise ValueError(msg)

    annotation_range_lock = annotation_run.annotation_range_lock
    if annotation_range_lock is None:
        msg = "Can't retry annotation run with no annotation lock!"
        raise ValueError(msg)

    annotation_run.error_exception = None  # Clear so status won't be error
    annotation_run.task_id = None  # Allow celery jobs to get lock on it
    annotation_run.save()

    if upload_only:
        # Delete uploaded data, then hand back to the dispatcher: the run stays ANNOTATION_COMPLETED with
        # its annotated VCF present -> import lane -> db_workers (import_annotation_run). #1649
        tasks = [
            delete_annotation_run_uploaded_data.si(annotation_run.pk),
        ]
    else:
        # Delete old AnnotationRun then try again: a fresh CREATED run -> VEP lane (annotate_variants). #1649
        old_annotation_run = annotation_run
        annotation_run = AnnotationRun.objects.create(pipeline_type=old_annotation_run.pipeline_type)
        tasks = [
            delete_annotation_run.si(old_annotation_run.pk),
            assign_range_lock_to_annotation_run.si(annotation_run.pk, annotation_range_lock.pk),
        ]

    # #1649: retry no longer launches annotate_variants inline - it hands back to the single-authority
    # dispatcher (preserving lease/reclaim/attempt-cap), which launches the right lane by status.
    tasks.append(_dispatch_trigger_sig(annotation_range_lock.version_id))
    task = chain(tasks)
    task.apply_async()

    return annotation_run


def _unannotated_variants_to_vcf(genome_build: GenomeBuild, vcf_filename,
                                 annotation_range_lock, pipeline_type) -> int:
    logging.info("unannotated_variants_to_vcf()")
    if os.path.exists(vcf_filename):
        raise ValueError(f"Don't want to overwrite '{vcf_filename}' which already exists!")
    mk_path_for_file(vcf_filename)
    kwargs = {}
    if annotation_range_lock:
        kwargs["min_variant_id"] = annotation_range_lock.min_variant.pk
        kwargs["max_variant_id"] = annotation_range_lock.max_variant.pk

    annotation_version = annotation_range_lock.version.get_any_annotation_version()
    qs = get_variants_qs_for_annotation(annotation_version, pipeline_type=pipeline_type, **kwargs)

    if pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
        if settings.ANNOTATION_VEP_SV_MAX_SIZE:
            # VEP will skip variants above a certain size and fill up the logs with 'too long to annotate'
            # So just skip these. I don't think it makes much difference in memory usage
            q_not_too_long = Q(svlen__isnull=True) | Q(abs_svlen__lte=settings.ANNOTATION_VEP_SV_MAX_SIZE)
            qs = qs.annotate(abs_svlen=Abs("svlen")).filter(q_not_too_long)
    return write_qs_to_vcf(vcf_filename, genome_build, qs)


def write_qs_to_vcf(vcf_filename, genome_build, qs, info_dict=VARIANT_GRID_INFO_DICT, use_accession=False) -> int:
    # We had an issue with writing accessions in VEP, so use chrom names and the default VEP fasta instead
    # @see https://github.com/Ensembl/ensembl-vep/issues/1635
    qs = qs.order_by("locus__contig__genomebuildcontig__order", "locus__position")
    if use_accession:
        chrom_key = "locus__contig__refseq_accession"
    else:
        chrom_key = "locus__contig__name"

    sorted_values = qs.values("id", chrom_key, "locus__position",
                              "locus__ref__seq", "alt__seq", "end", "svlen")

    # External dumps are written .vcf.gz (see AnnotationRun.get_dump_filename) - gzip on the way out
    if vcf_filename.endswith(".gz"):
        f = gzip.open(vcf_filename, "wt", compresslevel=6)
    else:
        f = open(vcf_filename, "wt")
    with f:
        return write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict=info_dict,
                                                      use_accession=use_accession)
