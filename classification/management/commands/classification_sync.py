from django.core.management import BaseCommand

from sync.models import SyncDestination
from sync.sync_run import run_sync


class Command(BaseCommand):
    """
    Performs a sync to other VariantGrid instances
    """
    category = "ops"

    def add_arguments(self, parser):
        parser.add_argument('--destination', required=True)

    def handle(self, *args, **options):
        destination_str = options["destination"]
        sd: SyncDestination = SyncDestination.objects.get(name=destination_str)
        run_sync(sd, full_sync=True)
