from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.upgrader import record_attempt


class Command(BaseCommand):
    """ Record an attempt at a manual task by hand (the upgrader records its own) """
    category = "ops"

    def add_arguments(self, parser: CommandParser):
        parser.add_argument('--id', required=True, help="Command of the id that has been completed")
        parser.add_argument('--failed', action='store_true', help="If attempted but failed")
        parser.add_argument('--note', help="Optional note")
        parser.add_argument('--ver', help="Version of the code this was run against")

    def handle(self, *args, **options):
        record_attempt(options["id"], success=not options["failed"], note=options.get("note"),
                       version=options.get("ver"))
        print("Attempt added")
