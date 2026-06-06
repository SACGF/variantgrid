"""
Celery task that performs the long-running archive of a VCF's underlying data.

The precondition (uploaded source file present) is validated synchronously in the
web request; the slow work (zygosity count walk + dropping partition data) runs
here so the request doesn't time out.

@see snpdb.archive
"""

import logging

import celery
from django.contrib.auth.models import User

from snpdb.archive import archive_vcf
from snpdb.models import VCF


@celery.shared_task(queue="db_workers")
def archive_vcf_task(vcf_id: int, user_id: int, reason: str = "", force: bool = False):
    vcf = VCF.objects.get(pk=vcf_id)
    user = User.objects.get(pk=user_id)
    logging.info("archive_vcf_task: archiving VCF %s (force=%s)", vcf_id, force)
    try:
        archive_vcf(vcf, user, reason=reason, force=force)
    except Exception:
        # archive_vcf rolled back its transaction; clear the in-progress marker (set
        # before queueing) so the VCF can be archived again.
        VCF.objects.filter(pk=vcf_id).update(data_archive_started_date=None)
        raise
