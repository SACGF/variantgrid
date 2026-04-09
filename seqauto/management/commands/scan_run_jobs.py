from django.core.management.base import BaseCommand

from seqauto.models import SeqAutoRun
from seqauto.tasks.scan_run_jobs import process_seq_auto_run


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--process-types', action='append')
        parser.add_argument('--launch-types', action='append')
        parser.add_argument('--run-launch-script', action='store_true')
        parser.add_argument('--reuse-prev-scan-id', type=int, required=False,
                            help="Reused previous scanned files from run ID")

    def handle(self, *args, **options):
        only_process_file_types = options.get("process_types")
        only_launch_file_types = options.get("launch_types")
        run_launch_script = options.get("run_launch_script", False)
        reuse_prev_scan_id = options.get("reuse_prev_scan_id")

        seqauto_run = SeqAutoRun.objects.create()
        process_seq_auto_run(seq_auto_run_id=seqauto_run.pk,  # @UndefinedVariable
                             only_process_file_types=only_process_file_types,
                             only_launch_file_types=only_launch_file_types,
                             run_launch_script=run_launch_script,
                             reuse_prev_scan_id=reuse_prev_scan_id)
