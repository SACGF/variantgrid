import logging

import celery
from django.conf import settings
from django.core.cache import cache

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.models import AnnotationRun, VariantAnnotationPipelineType
from annotation.models.models import AnnotationVersion, AnnotationRangeLock
from annotation.tasks.annotate_variants import annotate_variants
from library.log_utils import log_traceback
from snpdb.models import GenomeBuild, ImportStatus, Sample, VCF, Variant


@celery.shared_task(queue='scheduling_single_worker')
def annotation_scheduler(active=True):
    """ This is run on scheduling_single_worker queue to avoid race conditions """
    LOCK_EXPIRE = 60 * 5  # 5 minutes
    lock_id = "annotation-scheduler-lock"

    # cache.add fails if the key already exists
    acquire_lock = lambda: cache.add(lock_id, "true", LOCK_EXPIRE)
    release_lock = lambda: cache.delete(lock_id)

    try:
        if acquire_lock():
            try:
                logging.info("Got the lock for annotation scheduler")
                for genome_build in GenomeBuild.builds_with_annotation():
                    annotation_version = AnnotationVersion.latest(genome_build, active=active)
                    variant_annotation_version = annotation_version.variant_annotation_version
                    while True:
                        range_lock = _handle_variant_annotation_version(variant_annotation_version)
                        if range_lock is None:
                            break
            finally:
                logging.info("Releasing lock")
                release_lock()
        else:
            logging.info("Someone else has %s", lock_id)
    except:
        log_traceback()


def _handle_range_lock(range_lock, pipeline_type=None):
    pipeline_types = []
    if pipeline_type is not None:
        pipeline_types.append(pipeline_type)
    else:
        pipeline_types = list(VariantAnnotationPipelineType)

    for pipeline_type in pipeline_types:
        annotation_run, created = AnnotationRun.objects.get_or_create(annotation_range_lock=range_lock,
                                                                      pipeline_type=pipeline_type)
        if created:
            annotate_variants.apply_async((annotation_run.pk,))  # @UndefinedVariable


def _handle_variant_annotation_version(variant_annotation_version):
    # If we crash in the wrong place, we may end up with unassigned range locks (need 1 of each type)
    arl_qs = AnnotationRangeLock.objects.filter(version=variant_annotation_version)
    for pipeline_type in VariantAnnotationPipelineType:
        # Look for missing any AnnotationRun for this lock
        for range_lock in arl_qs.exclude(annotationrun__pipeline_type=pipeline_type):
            logging.warning("Assigned orphaned annotation range lock: %s, run type: %s",
                            range_lock, VariantAnnotationPipelineType(pipeline_type).label)
            _handle_range_lock(range_lock, pipeline_type)

    range_lock, unannotated_count = get_annotation_range_lock_and_unannotated_count(variant_annotation_version,
                                                                                    settings.ANNOTATION_VEP_BATCH_MIN,
                                                                                    settings.ANNOTATION_VEP_BATCH_MAX)
    if range_lock is not None:
        range_lock.save()
        _handle_range_lock(range_lock)
    else:
        if unannotated_count:
            logging.warning("Unannotated variants but couldn't get a lock.")
            logging.warning("If you are sure no annotation requests are running you can reset things with 'reset_annotation_states'")
        else:
            logging.info("No unannotated variants left!")

            waiting_vcfs = VCF.objects.filter(import_status=ImportStatus.IMPORTING)
            waiting_samples = Sample.objects.filter(import_status=ImportStatus.IMPORTING)
            if waiting_vcfs or waiting_samples:
                logging.warning("Waiting vcfs or samples - if you are sure no annotation requests are running you can set vcfs to be available with 'annotation_set_all_complete' ")
            else:
                logging.info("No waiting vcfs or samples (this is good!)")

    return range_lock


def subdivide_annotation_range_lock(arl: AnnotationRangeLock,
                                    minimum_size=AnnotationRangeLock.MIN_SIZE_FOR_SUBDIVISION) -> AnnotationRangeLock:
    """ Sometimes we have an annotation run crash or timeout etc, this splits it in half so we can try again
        With a smaller run

        Range lock passed in becomes bottom half of original range
        Returns new range lock for top half of original range
    """

    size = int(arl.max_variant_id) - int(arl.min_variant_id)
    if size < minimum_size:
        raise ValueError(f"Cannot subdivide {arl} below minimum size of {minimum_size}")

    logging.info("Subdividing %s", arl)

    # Delete all existing annotation runs
    res = arl.annotationrun_set.all().delete()
    logging.info("Deleted attached annotation runs: %s", res)

    half_size = size / 2
    halfway_point = int(arl.min_variant_id) + half_size
    first_at_or_above_halfway = Variant.objects.filter(pk__gte=halfway_point).order_by("pk").first()
    first_below_halfway = Variant.objects.filter(pk__gte=arl.min_variant_id, pk__lt=halfway_point).order_by("pk").last()

    av = arl.version.get_any_annotation_version()
    unannotated_qs = get_variants_qs_for_annotation(av, min_variant_id=first_at_or_above_halfway.pk,
                                                    max_variant_id=arl.max_variant.pk)

    top_half_count = unannotated_qs.count()
    new_lock = AnnotationRangeLock.objects.create(version=arl.version,
                                                  min_variant=first_at_or_above_halfway,
                                                  max_variant=arl.max_variant,
                                                  count=top_half_count)

    arl.max_variant = first_below_halfway
    arl.count -= top_half_count
    arl.save()

    logging.info("Shifted old lock down: %s", arl)
    logging.info("New lock created: %s", new_lock)
    return new_lock
