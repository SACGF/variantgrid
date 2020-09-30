from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.models import ManualMigrationTask, ManualMigrationRequired


class Command(BaseCommand):
    """
    To be called via the migrator.py program
    """

    def add_arguments(self, parser: CommandParser):
        parser.add_argument('--category', required=True, choices=['manage', 'pip', 'admin', 'other'],
                            help="The kind of task, manage for manage.py, pip for requirement, admin for a task to be run from an admin screen, otherwise other")
        parser.add_argument('--line', required=True, help="The task that has to take place, g.g. the text after python3.8 manage.py")
        parser.add_argument('--note', required=False, help="Optional human friendly note")

    def handle(self, *args, **options):
        category = options["category"]
        line = options["line"]
        note = options.get("note")

        task_id = f"{category}*{line}"
        mmt, _ = ManualMigrationTask.objects.get_or_create(id=task_id)

        ManualMigrationRequired.objects.create(
            task=mmt,
            note=note
        )

        print("Requirement added")
