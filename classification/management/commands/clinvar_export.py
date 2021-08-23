from django.core.management import BaseCommand

from classification.models.clinvar_export_prepare import ClinvarAlleleExportPrepare
from snpdb.models import Allele
from sync.models import SyncDestination


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--prepare',  action='store_true', default=False)

    def handle(self, *args, **options):
        if options["prepare"]:
            for count, allele in enumerate(Allele.objects.all()):
                report = ClinvarAlleleExportPrepare(allele).update_export_records()
                if count % 100 == 0:
                    print(f"Processing Allele no {count}")
                # print(report)
            print(f"Completed {count}")
