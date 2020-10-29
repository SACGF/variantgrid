from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from classification.models import ConditionAlias


class Command(BaseCommand):
    """
    Updates ConditionAlias to match the state of the classifications
    """

    def handle(self, *args, **options):
        ConditionAlias.sync_aliases()