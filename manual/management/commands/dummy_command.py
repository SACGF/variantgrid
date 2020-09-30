from django.core.management import CommandParser
from django.core.management.base import BaseCommand


class Command(BaseCommand):

    def add_arguments(self, parser: CommandParser):
        pass

    def handle(self, *args, **options):
        print("Dummy command passed")
