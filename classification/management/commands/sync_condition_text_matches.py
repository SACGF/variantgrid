from django.core.management.base import BaseCommand
from classification.models.condition_text_matching import ConditionTextMatch


class Command(BaseCommand):
    """
    Updates Condition Text Matches to match the state of the classifications
    """

    def handle(self, *args, **options):
        print("Syncing")
        ConditionTextMatch.sync_all()
        print("Complete")