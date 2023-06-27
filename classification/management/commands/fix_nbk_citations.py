from django.core.management import BaseCommand
from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--delete', action='store_true', default=False)

    def handle(self, *args, **options):
        for classification in Classification.objects.iterator():
            has_nbk = False
            for ref in classification.db_refs:
                if ref.get("db") == "HTTPS" and "www.ncbi.nlm.nih.gov/books" in ref.get("id"):
                    has_nbk = True
                    break

            if has_nbk:
                print(f"Revalidating classification {classification.pk}")
                classification.revalidate(user=admin_bot())
