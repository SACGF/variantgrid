import itertools

from django.core.management import BaseCommand

from classification.classification_import import variant_matching_dry_run
from classification.models import ImportedAlleleInfo


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        print("Triggering Allele Info validation, this may take a while")
        variant_matching_dry_run(ImportedAlleleInfo.objects.all())
        print("Complete")
