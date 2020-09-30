#!/usr/bin/env python3

from django.core.management.base import BaseCommand

from library.log_utils import console_logger
from seqauto.models import JobScript, JobScriptStatus


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--job_script_id', required=True)
        parser.add_argument('--return_code', type=int, required=True)

    def handle(self, *args, **options):
        job_script_pk = options["job_script_id"]
        return_code = options["return_code"]

        logger = console_logger()
        logger.info("Setting job %s to %d", job_script_pk, return_code)
        job_script = JobScript.objects.get(pk=job_script_pk)

        if job_script.job_status == JobScriptStatus.FINISHED:
            msg = f"JobScript {job_script_pk} already set to FINISHED!"
            logger.warning(msg)
        #    raise CommandError(msg)

        job_script.return_code = return_code
        job_script.job_complete()

        # TODO: What about dependent jobs which will be cancelled?
