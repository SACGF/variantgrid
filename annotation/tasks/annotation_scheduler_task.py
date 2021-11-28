import logging

import celery
from django.conf import settings
from django.core.cache import cache

from annotation.annotation_versions import get_annotation_range_lock_and_unannotated_count
from annotation.models import AnnotationRun
from annotation.models.models import AnnotationVersion
from annotation.tasks.annotate_variants import annotate_variants
from library.log_utils import log_traceback
from snpdb.models import GenomeBuild, ImportStatus, Sample, VCF


@celery.shared_task
def annotation_scheduler():
    """ This is run on scheduling_single_worker queue to avoid race conditions """
    LOCK_EXPIRE = 60 * 5  # 5 minutes
    lock_id = "annotation-scheduler-lock"

    # cache.add fails if if the key already exists
    acquire_lock = lambda: cache.add(lock_id, "true", LOCK_EXPIRE)
    release_lock = lambda: cache.delete(lock_id)

    try:
        if acquire_lock():
            try:
                logging.info("Got the lock for annotation scheduler")
                for genome_build in GenomeBuild.builds_with_annotation():
                    annotation_version = AnnotationVersion.latest(genome_build)
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


def _handle_variant_annotation_version(variant_annotation_version):
    range_lock, unannotated_count = get_annotation_range_lock_and_unannotated_count(variant_annotation_version,
                                                                                    settings.ANNOTATION_VEP_BATCH_MIN,
                                                                                    settings.ANNOTATION_VEP_BATCH_MAX)
    if range_lock is not None:
        range_lock.save()
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=range_lock)
        annotate_variants.apply_async((annotation_run.pk,))  # @UndefinedVariable
        logging.info("Scheduled annotate_variants job")
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
