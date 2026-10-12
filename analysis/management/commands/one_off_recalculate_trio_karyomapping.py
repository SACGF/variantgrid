"""
Recount every trio's GenomeKaryomappingCounts. Joint-called trios store a reference call as HOM_REF, which the
karyotype bins used to miss (they only matched MISSING), leaving those counts near zero (#1754).
One celery task per trio, as each reads the whole cohort.
"""
from django.core.management import BaseCommand

from analysis.tasks.karyomapping_tasks import create_genome_karyomapping_for_trio
from snpdb.models import Trio


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        num_queued = 0
        for trio in Trio.objects.select_related("cohort").order_by("pk"):
            if trio.cohort.data_archived:
                continue
            create_genome_karyomapping_for_trio.si(trio.pk).apply_async()
            num_queued += 1
        self.stdout.write(f"Queued karyomapping for {num_queued} trios")
