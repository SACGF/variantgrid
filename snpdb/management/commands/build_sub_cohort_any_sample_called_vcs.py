from django.core.management.base import BaseCommand

from snpdb.models import Cohort, ImportStatus
from snpdb.tasks.sub_cohort_tasks import build_sub_cohort_any_sample_called_vc_task


class Command(BaseCommand):
    """ Backfill the pre-computed any-sample-called VariantCollection for existing sub-cohorts (issue #1551).

        Runs the build synchronously so the operator can monitor progress, and skips sub-cohorts that
        already have a current collection. """

    def add_arguments(self, parser):
        parser.add_argument('--force', action='store_true',
                            help="Rebuild even if a current collection already exists")

    def handle(self, *args, **options):
        force = options["force"]
        qs = Cohort.objects.filter(parent_cohort__isnull=False,
                                   import_status=ImportStatus.SUCCESS).order_by("pk")
        total = qs.count()
        print(f"Building any-sample-called VariantCollections for {total} sub-cohorts...")
        built = 0
        skipped = 0
        for i, cohort in enumerate(qs):
            if not force and cohort.get_any_sample_called_variant_collection() is not None:
                skipped += 1
                continue
            print(f"{100 * i / total:.1f}% - building for sub-cohort {cohort.pk}: {cohort}")
            build_sub_cohort_any_sample_called_vc_task(cohort.pk)  # synchronous
            built += 1
        print(f"Done. built={built}, skipped(already current)={skipped}")
