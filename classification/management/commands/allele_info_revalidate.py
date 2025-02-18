from django.core.management import BaseCommand

from classification.classification_import import variant_matching_dry_run
from classification.models import ImportedAlleleInfo


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument("--dirty", action="store_true", help="Only validate Allele Infos already marked as dirty")

    def handle(self, *args, **options):
        qs = ImportedAlleleInfo.objects.all()
        if options["dirty"]:
            qs = qs.filter(dirty_message__isnull=False)
        print("Triggering Allele Info validation, this may take a while")
        variant_matching_dry_run(qs)
        print("Complete")
