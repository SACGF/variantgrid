import logging

import celery
from django.db.models.aggregates import Count

from library.guardian_utils import check_can_write
from library.log_utils import log_traceback
from snpdb.models import VCF, ImportStatus, Sample
from snpdb.variant_zygosity_count import update_all_variant_zygosity_counts_for_vcf, \
    update_all_variant_zygosity_counts_for_sample


def soft_delete_vcfs(user, *vcf_ids):
    vcf_ids_to_delete = []

    for vcf_id in vcf_ids:
        vcf = VCF.get_for_user(user, vcf_id)
        check_can_write(vcf, user)
        vcf_ids_to_delete.append(vcf_id)

    logging.info("Going to delete vcfs: %s", ','.join(map(str, vcf_ids_to_delete)))

    # Soft-Delete vcfs and samples - hide from grids/lists etc and delete async
    vcfs_marked_for_deletion = VCF.objects.filter(pk__in=vcf_ids_to_delete).update(import_status=ImportStatus.MARKED_FOR_DELETION)
    samples_marked_for_deletion = Sample.soft_delete_samples_with_deleted_vcfs()
    logging.info("Marked %d vcfs, %d samples for deletion", vcfs_marked_for_deletion, samples_marked_for_deletion)
    remove_soft_deleted_vcfs_task.apply_async(countdown=1)  # To make sure that vcfs have been set to deleted


@celery.shared_task(ignore_result=True)
def remove_soft_deleted_vcfs_task():
    """ This is a clean up job so only want to run 1 copy - ie on schedule_single_worker queue

        If this becomes a bottleneck, we could run update_variant_zygosity_count tasks with deletion task afterwards.
    """

    try:

        # Note: It's better to do this one vcf at a time as it appears to blow out commit logs
        deleted_vcfs = 0

        # If anything was left in DELETING state, since this is only process - can just try again.
        proj_potentials = VCF.objects.filter(import_status__in=ImportStatus.DELETION_STATES).values_list('pk', flat=True)
        for pk in proj_potentials:
            qs = VCF.objects.filter(pk=pk)
            qs.filter(import_status=ImportStatus.MARKED_FOR_DELETION).update(import_status=ImportStatus.DELETING)
            vcf = qs.get()
            logging.info("Deleting vcfs %s", vcf)
            try:
                update_all_variant_zygosity_counts_for_vcf(vcf, '-')
            except:  # Already logged etc
                pass

            vcf.delete()
            deleted_vcfs += 1

        logging.info("Deleted %d vcfs", deleted_vcfs)

        deleted_samples = 0

        # Only handle samples where vcfs are not being deleted (full vcf deletions happen above)
        # May have been scheduled
        single_sample_deletions_qs = Sample.objects.exclude(vcf__import_status__in=ImportStatus.DELETION_STATES)
        single_sample_deletions_qs = single_sample_deletions_qs.filter(import_status=ImportStatus.MARKED_FOR_DELETION)
        ss_potentials = single_sample_deletions_qs.values_list('pk', flat=True)
        logging.info("remove_soft_deleted_vcfs_task - Found %d samples marked for deletion", ss_potentials.count())

        for pk in ss_potentials:
            qs = Sample.objects.filter(pk=pk)
            qs.filter(import_status=ImportStatus.MARKED_FOR_DELETION).update(import_status=ImportStatus.DELETING)
            sample = qs.get()
            logging.info("remove_soft_deleted_vcfs_task - deleting sample %s", sample)
            try:
                update_all_variant_zygosity_counts_for_sample(sample, '-')
            except:  # Already logged etc
                pass
            sample.delete()
            deleted_samples += 1

        vcfs_with_no_samples_qs = VCF.objects.all().annotate(sample_count=Count("sample")).filter(sample_count=0)
        deleted_empty_vcfs = vcfs_with_no_samples_qs.count()
        if deleted_empty_vcfs:
            logging.info("Deleting %d vcfs with no samples (no count adjustments should be needed)", deleted_empty_vcfs)
            vcfs_with_no_samples_qs.delete()

        logging.info("Deleted %d samples (then %d empty vcfs)", deleted_samples, deleted_empty_vcfs)
    except:
        log_traceback()
