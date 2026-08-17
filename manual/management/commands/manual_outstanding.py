import json

from django.core.management import CommandParser, get_commands
from django.core.management.base import BaseCommand

from manual.gates import blocked_by
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
        known_commands = set(get_commands())
        outstanding_tasks = ManualMigrationOutstanding.outstanding_tasks()
        task_list = []
        for outstanding_task in outstanding_tasks:
            task_json = outstanding_task.to_json()
            is_manage = task_json["category"] == "manage"
            # command_exists only applies to 'manage' tasks (None for 'other' human steps, which have no
            # command). A manage task whose command was deleted is an obsolete/orphaned row - never
            # auto-run it (it would crash), it's reported for review instead.
            command_exists = task_json["line"].split()[0] in known_commands if is_manage else None
            task_json["command_exists"] = command_exists
            blocking = blocked_by(task_json.get("requires"))
            task_json["blocked_by"] = blocking
            # Only 'manage' tasks are auto-runnable, and only once nothing blocks them and the command exists.
            task_json["runnable"] = bool(is_manage and command_exists and not blocking)
            task_list.append(task_json)
        envelope = {"tasks": task_list,
                    "ROLLBAR_ACCESS_TOKEN": rollbar_token}
        print(json.dumps(envelope))
