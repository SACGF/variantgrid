from django.core.management.base import BaseCommand

from classification.models import ConditionText, sync_all_condition_resolutions_to_classifications
from classification.models.condition_text_matching import ConditionTextMatch


class Command(BaseCommand):
    """
    Updates Condition Text Matches to match the state of the classifications
    """

    def add_arguments(self, parser):
        parser.add_argument('--reset', action='store_true', default=False)
        parser.add_argument('--clear', action='store_true', default=False)
        parser.add_argument('--classifications', action='store_true', default=False)

    def handle(self, *args, **options):
        if options["classifications"]:
            print("Updating classifications")
            sync_all_condition_resolutions_to_classifications()
            print("Complete")
            return

        if options["reset"]:
            print("Deleting old records")
            ConditionText.objects.all().delete()
        if options["clear"]:
            print("Clearing all existing values")
            ct: ConditionText
            for ct in ConditionText.objects.all():
                ct.clear()

        print("Syncing")
        ConditionTextMatch.sync_all()
        print("Complete")