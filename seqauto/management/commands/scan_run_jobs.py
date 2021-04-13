from django.core.management.base import BaseCommand
from seqauto.tasks.scan_run_jobs import scan_run_jobs


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--process-types', action='append')
        parser.add_argument('--launch-types', action='append')
        parser.add_argument('--run-launch-script', action='store_true')

    def handle(self, *args, **options):
        only_process_file_types = options.get("process_types")
        only_launch_file_types = options.get("launch_types")
        run_launch_script = options.get("run_launch_script", False)

        scan_run_jobs(only_process_file_types=only_process_file_types,  # @UndefinedVariable
                      only_launch_file_types=only_launch_file_types,
                      run_launch_script=run_launch_script)
