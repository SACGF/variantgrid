import json

from django.core.management.base import BaseCommand

from manual.upgrader import outstanding_tasks_json


class Command(BaseCommand):
    """ Outstanding manual tasks as JSON, with whether each is runnable / blocked / obsolete """
    category = "ops"

    def handle(self, *args, **options):
        print(json.dumps({"tasks": outstanding_tasks_json()}))
