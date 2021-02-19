from django.core.management.base import BaseCommand

from classification.models import ConditionText
from classification.models.condition_text_matching import ConditionTextMatch


class Command(BaseCommand):
    """
    Updates Condition Text Matches to match the state of the classifications
    """

    def add_arguments(self, parser):
        parser.add_argument('--reset', action='store_true', default=False)

    def handle(self, *args, **options):
        if options["reset"]:
            print("Deleting old records")
            ConditionText.objects.all().delete()

        force = options["force"]

        print("Syncing")
        ConditionTextMatch.sync_all()
        print("Complete")