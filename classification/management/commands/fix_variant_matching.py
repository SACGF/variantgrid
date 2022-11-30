import csv
import time
from time import sleep
from typing import List, Dict, Optional, Set

from django.core.management import BaseCommand
from django.db.models import Q

from classification.classification_import import process_classification_import
from classification.enums import SpecialEKeys
from classification.models import Classification, ClassificationImport
from classification.models.classification_import_run import ClassificationImportRun
from library.guardian_utils import admin_bot
from snpdb.models import ImportSource


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', default=False, help='Attempt to rematch every single classification')
        parser.add_argument('--file', type=str, help='If provided expects a csv where the first column is a combination of genome build and imported c.hgvs it expects to import')
        parser.add_argument('--missing', action='store_true', default=False, help='Attempt to rematch only classifications not linked to a variant - one at a time')
        parser.add_argument('--extra', action='store_true', default=False, help='Populate the variant_info of a classification')


    def report_unmatched(self):
        print(f"Unmatched count = {Classification.objects.filter(variant__isnull=True).count()}")

    def handle(self, *args, **options):
        mode_all = options.get('all')
        mode_missing = options.get('missing')
        mode_extra = options.get('extra')
        mode_file = options.get('file')

        if mode_all and mode_missing:
            raise ValueError("all and missing are mutually exclusive parameters")

        if mode_file and (mode_all or mode_missing or mode_extra):
            raise ValueError("file is an exclusive parameter")

        if not any([mode_all, mode_missing, mode_file, mode_extra]):
            raise ValueError("Need one of all, file, missing, mode_extra")

        row_count = 0
        batch_size = 50
        import_keys: Optional[Set[str]] = None

        if mode_extra:
            self.handle_extra()
            return
        if mode_missing:
            batch_size = 1
        elif mode_file:
            header_row = True
            import_keys = set()
            with open(mode_file, 'r') as file_handle:
                for row in csv.reader(file_handle):
                    if header_row:
                        if not row[0] == "Import Key":
                            raise ValueError("Expected first column of file to be - 'Import Key'")
                        header_row = False
                    else:
                        import_keys.add(row[0].strip())
            print(f"Only running over import keys in file - {len(import_keys)}")

        self.report_unmatched()

        qs = Classification.objects.all().order_by('evidence__genome_build__value')
        if mode_missing:
            qs = qs.filter(variant__isnull=True)

        user = admin_bot()

        # setup a temporary import so discordance notifications are not sent out
        try:
            batch: List[Classification] = []
            for c in qs:
                if import_keys is not None:
                    import_key = f"{c.get(SpecialEKeys.GENOME_BUILD) or ''}#{c.get(SpecialEKeys.C_HGVS) or ''}"
                    if import_key not in import_keys:
                        continue

                c.revalidate(user=user)
                batch.append(c)
                row_count += 1

                if len(batch) >= batch_size:
                    self.handle_batch(batch)
                    batch = []
                    print(f"Handled {row_count}")

            self.handle_batch(batch)
            print(f"Handled {row_count}")
            self.report_unmatched()
            sleep(10)  # give time for variant matching to complete
        finally:
            ClassificationImportRun.record_classification_import("variant_rematching", 0, is_complete=True)

    def sleep_for_delay(self):
        time.sleep(10)

    def handle_extra(self):
        for i, c in enumerate(Classification.objects.all()):
            if i % 100 == 0:
                print(f"Processed {i} classifications")
            c.update_allele_info()
            c.save(update_modified=False)
        print(f"Finished {i} classifications")


    def handle_batch(self, batch: List[Classification]):
        ClassificationImportRun.record_classification_import("variant_rematching", len(batch))
        user = admin_bot()
        if batch:
            imports_by_genome: Dict[int, ClassificationImport] = {}
            for vc in batch:
                try:
                    genome_build = vc.get_genome_build()
                    if genome_build.pk not in imports_by_genome:
                        imports_by_genome[genome_build.pk] = ClassificationImport.objects.create(user=user,
                                                                                                 genome_build=genome_build)
                    vc_import = imports_by_genome[genome_build.pk]
                    vc.set_variant(variant=None, message='Admin has re-triggered variant matching')
                    vc.classification_import = vc_import
                    vc.save()
                except ValueError as ve:
                    print(f"Couldn't revalidate {vc.id} due to bad genome build {ve}")

            for vc_import in imports_by_genome.values():
                process_classification_import(vc_import, ImportSource.API)

            self.sleep_for_delay()
