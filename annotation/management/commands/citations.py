from django.core.management import BaseCommand

from annotation.models import Citation, CitationFetchRequest
from library.utils import batch_iterator


class Command(BaseCommand):

    def add_arguments(self, parser):
        #parser.add_argument('--all', action='store_true', default=False, help='Attempt to rematch every single classification')
        pass

    def handle(self, *args, **options):
        unloaded_citations = Citation.objects.filter(last_loaded__isnull=True)
        print(f"{unloaded_citations.count()} unloaded citations")
        loaded = 0
        for batch in batch_iterator(unloaded_citations, batch_size=20):
            CitationFetchRequest.fetch_all_now(batch)
            loaded += len(batch)
            print(f"Loaded {loaded} citations")
