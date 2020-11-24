from django.conf import settings
from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from library.git import Git
from manual.models.deployment_models import Deployment


class Command(BaseCommand):

    def add_arguments(self, parser: CommandParser):
        pass

    def handle(self, *args, **options):
        git_hash = Git(settings.BASE_DIR).hash
        Deployment.objects.create(git_hash=git_hash)
        print("Deployment recorded")
