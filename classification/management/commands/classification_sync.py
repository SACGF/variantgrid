from django.core.management import BaseCommand

from sync.models import SyncDestination


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--destination', required=True)

    def handle(self, *args, **options):
        destination_str = options["destination"]
        sd: SyncDestination = SyncDestination.objects.get(name=destination_str)
        sd.run(full_sync=True)
