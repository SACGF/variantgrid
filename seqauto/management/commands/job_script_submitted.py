#!/usr/bin/env python3

from django.core.management.base import BaseCommand

from seqauto.models import JobScript, JobScriptStatus


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--job_script_id', type=int, required=True)
        parser.add_argument('--job_id')

    def handle(self, *args, **options):
        job_script_id = options["job_script_id"]
        job_id = options.get("job_id")

        job_script = JobScript.objects.get(id=job_script_id)
        job_script.job_status = JobScriptStatus.SUBMITTED
        if job_id:
            job_script.job_id = job_id
        job_script.save()
