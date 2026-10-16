"""
Backfill ResolvedVariantInfo.p_hgvs / p_hgvs_source / p_hgvs_error (#1383). Only those three columns are written,
so it is cheap and safe to rerun - eg after a new VEP version lands or transcript sequences have been cached.

Rows are split into pk ranges of --batch-size and queued as celery tasks on db_workers, so batches run in parallel
(--sync runs them in this process instead).

The default pass never fetches transcript sequences; --retrieve-transcripts is a second, slower pass that fetches
the missing ones from NCBI / Ensembl for rows the first pass couldn't translate (use with --only-missing).
"""
from django.core.management.base import BaseCommand
from django.db.models import Count

from classification.models.classification_variant_info_models import ResolvedVariantInfo
from classification.tasks.resolved_variant_info_p_hgvs_task import (
    pk_ranges,
    recalc_p_hgvs_for_pk_range,
    recalc_p_hgvs_for_pk_range_task,
    resolved_variant_infos_for_p_hgvs,
)


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', help="Only this build (default all)")
        parser.add_argument('--only-missing', action='store_true', help="Only rows with no p_hgvs")
        parser.add_argument('--retrieve-transcripts', action='store_true',
                            help="Fetch missing transcript sequences (slow, network)")
        parser.add_argument('--batch-size', type=int, default=1000)
        parser.add_argument("--sync", action="store_true", help="Run in this process rather than queuing celery tasks")

    def handle(self, *args, **options):
        batch_kwargs = {
            "genome_build_name": options["genome_build"],
            "only_missing": options["only_missing"],
            "retrieve_transcripts": options["retrieve_transcripts"],
        }
        rvi_qs = resolved_variant_infos_for_p_hgvs(options["genome_build"], options["only_missing"])
        ranges = list(pk_ranges(rvi_qs, options["batch_size"]))

        if options["sync"]:
            num_processed = sum(recalc_p_hgvs_for_pk_range(min_pk, max_pk, **batch_kwargs)
                                for min_pk, max_pk in ranges)
            self.stdout.write(f"Processed {num_processed} ResolvedVariantInfo rows")
            self._print_summary(resolved_variant_infos_for_p_hgvs(options["genome_build"]))
        else:
            for min_pk, max_pk in ranges:
                recalc_p_hgvs_for_pk_range_task.si(min_pk, max_pk, **batch_kwargs).apply_async()
            self.stdout.write(f"Queued {len(ranges)} batches on db_workers")

    def _print_summary(self, rvi_qs):
        source_labels = dict(ResolvedVariantInfo._meta.get_field("p_hgvs_source").choices)
        self.stdout.write("Rows per p.HGVS source:")
        for row in rvi_qs.values("p_hgvs_source").annotate(count=Count("pk")).order_by("-count"):
            self.stdout.write(f"  {source_labels.get(row['p_hgvs_source'], 'None')}: {row['count']}")

        error_qs = rvi_qs.filter(p_hgvs__isnull=True, p_hgvs_error__isnull=False)
        self.stdout.write(f"Rows with an error: {error_qs.count()}")
        for row in error_qs.values("p_hgvs_error").annotate(count=Count("pk")).order_by("-count")[:10]:
            self.stdout.write(f"  {row['count']}: {row['p_hgvs_error']}")
