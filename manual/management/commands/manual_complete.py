from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.models import ManualMigrationTask, ManualMigrationAttempt


class Command(BaseCommand):
    """
    To be called via the migrator.py program
    """

    def add_arguments(self, parser: CommandParser):
        parser.add_argument('--id', required=True, help="Command of the id that has been completed")
        parser.add_argument('--failed', action='store_true', help="If attempted but failed")
        parser.add_argument('--note', help="Optional note")
        parser.add_argument('--ver', help="Version of the code this was run against")

    def handle(self, *args, **options):
        task_id = options["id"]
        success = not options.get("failed")
        note = options.get("note")
        version = options.get("ver")

        mmt, _ = ManualMigrationTask.objects.get_or_create(id=task_id)

        ManualMigrationAttempt.objects.create(
            task=mmt,
            note=note,
            source_version=version,
            requires_retry=not success
        )

        print("Attempt added")
