from django.core.management import BaseCommand
from snpdb.models import Allele

class Command(BaseCommand):

    def handle(self, *args, **options):
        update_count = 0

        for allele in Allele.objects.all():
            allele.validate()

            update_count += 1
            if update_count % 100 == 0:
                print(f"Completed {update_count}")
        print(f"Bulk Validation of Allele - completed")
