from django.conf import settings
from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from library.utils import get_git_hash
from manual.models.deployment_models import Deployment


class Command(BaseCommand):

    def add_arguments(self, parser: CommandParser):
        pass

    def handle(self, *args, **options):
        git_hash = get_git_hash(settings.BASE_DIR)
        Deployment.objects.create(git_hash=git_hash)
        print("Deployment recorded")
