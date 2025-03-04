from django.conf import settings
from django.core.management import BaseCommand

from snpdb.models import Lab


class Command(BaseCommand):
    """
    Removes aspects from the lab that would cause Slack notifications or an overlap of upload s3 locations.
    Leaves Clinvar keys in place (since ClinVar uploads can be disabled by settings)
    Leaves emails in place (since emails can be disabled by settings)
    """

    def add_arguments(self, parser):
        parser.add_argument('--environment', type=str, required=False)

    def handle(self, *args, **options):
        if settings.SITE_NAME != options["environment"]:
            print(f"Confirm by providing --environment \"{settings.SITE_NAME}\"")
            return

        lab: Lab
        for lab in Lab.objects.iterator():
            if lab.name == "Admin Lab":
                continue
            if lab.slack_webhook:
                lab.slack_webhook = ""
                lab.save(update_fields=["slack_webhook"])
                print(f"Updating {lab} to remove slack webhook")
            if lab.upload_location and lab.upload_location.startswith("s3://") and not lab.upload_location.endswith("_test"):
                # TODO would be nice to verify if this location exists
                print(f"Updating {lab} from {lab.upload_location} to {lab.upload_location}_test")
                lab.upload_location = f"{lab.upload_location}_test"
                lab.save(update_fields=["upload_location"])