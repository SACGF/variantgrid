import json

from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.models import ManualMigrationOutstanding

from variantgrid.settings.components.secret_settings import get_secret


class Command(BaseCommand):
    """
    To be called via the migrator.py program
    """

    def add_arguments(self, parser: CommandParser):
        pass

    def handle(self, *args, **options):
        rollbar_token = get_secret("ROLLBAR.access_token")
        outstanding_tasks = ManualMigrationOutstanding.outstanding_tasks()
        task_list = []
        for outstanding_task in outstanding_tasks:
            task_list.append(outstanding_task.to_json())
        envelope = {"tasks": task_list,
                    "ROLLBAR_ACCESS_TOKEN": rollbar_token}
        print(json.dumps(envelope))
