import time
from typing import List, Dict, Optional

from django.core.management import BaseCommand

from classification.classification_import import process_classification_import
from classification.models import Classification, ClassificationImport
from library.guardian_utils import admin_bot
from snpdb.models import ImportSource


class Command(BaseCommand):

    def report_unmatched(self):
        print(f"Unmatched count = {Classification.objects.filter(variant__isnull=True).count()}")

    def handle(self, *args, **options):
        row_count = 0

        batch_size = 50
        mode: Optional[str] = None
        if options.get('all'):
            mode = 'all'
        elif options.get('missing'):
            mode = 'missing'
            batch_size = 1

        if not mode:
            print("Must confirm running over all records with --all or --missing")
            return

        self.report_unmatched()

        qs = Classification.objects.all().order_by('evidence__genome_build__value')
        if mode == 'missing':
            qs = qs.filter(variant__isnull=True)

        user = admin_bot()

        batch: List[Classification] = list()
        for c in qs:
            c.revalidate(user=user)
            batch.append(c)
            row_count += 1

            if len(batch) >= batch_size:
                self.handle_batch(batch)
                batch = list()
                print(f"Handled {row_count}")

        self.handle_batch(batch)
        print(f"Handled {row_count}")
        self.report_unmatched()

    def sleep_for_delay(self):
        time.sleep(10)

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', default=False)
        parser.add_argument('--missing', action='store_true', default=False)

    def handle_batch(self, batch: List[Classification]):
        user = admin_bot()
        if batch:
            imports_by_genome: Dict[int, ClassificationImport] = dict()
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
