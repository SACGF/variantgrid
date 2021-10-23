import time
from typing import List, Dict

from django.core.management import BaseCommand

from classification.models import Classification, ClassificationImport
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def sleep_for_delay(self):
        time.sleep(10)

    @property
    def batch_size(self) -> int:
        return 50

    def add_arguments(self, parser):
        parser.add_argument('--all', action='store_true', default=False)

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
            self.sleep_for_delay()

    def handle(self, *args, **options):
        row_count = 0
        apply_all = options.get('all')
        if not apply_all:
            print("Must confirm running over all records with --all")

        batch_size = self.batch_size
        user = admin_bot()

        batch: List[Classification] = list()
        for c in Classification.objects.all().order_by('evidence__genome_build__value'):
            c.revalidate(user=user)
            batch.append(c)
            row_count += 1

            if len(batch) >= batch_size:
                self.handle_batch(batch)
                batch = list()
                print(f"Handled {row_count}")

        self.handle_batch(batch)
        print(f"Handled {row_count}")
