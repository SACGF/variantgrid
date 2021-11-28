import json

from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.models import ManualMigrationOutstanding


class Command(BaseCommand):
    """
    To be called via the migrator.py program
    """

    def add_arguments(self, parser: CommandParser):
        pass

    def handle(self, *args, **options):
        outstanding_tasks = ManualMigrationOutstanding.outstanding_tasks()
        task_list = list()
        for outstanding_task in outstanding_tasks:
            task_list.append(outstanding_task.to_json())
        envelope = {"tasks": task_list}
        print(json.dumps(envelope))
