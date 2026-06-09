import logging

from django.core.management.base import BaseCommand
from django.db import connection

from annotation.models.models import VariantAnnotation, VariantAnnotationVersion


class Command(BaseCommand):
    help = ("Backfill VariantAnnotation.spliceai_max_ds. Idempotent. "
            "By default only processes ACTIVE (latest/current) versions per genome build - "
            "historical versions are archived or use the slower DamageNode fallback. "
            "Use --all to process every partition.")

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group()
        group.add_argument("--vav-id", type=int,
                           help="Only process this VariantAnnotationVersion")
        group.add_argument("--all", action="store_true",
                           help="Process all versions (default: only ACTIVE/latest per genome build)")
        parser.add_argument("--skip-index", action="store_true",
                            help="Skip CREATE INDEX CONCURRENTLY step")
        parser.add_argument("--chunk-size", type=int, default=10_000,
                            help="pk-range chunk size for the per-partition UPDATE loop")

    def handle(self, *args, **options):
        qs = VariantAnnotationVersion.objects.all().order_by("pk")
        if vav_id := options.get("vav_id"):
            qs = qs.filter(pk=vav_id)
        elif not options["all"]:
            qs = qs.filter(status=VariantAnnotationVersion.Status.ACTIVE)

        chunk_size = options["chunk_size"]
        for vav in qs:
            updated = VariantAnnotation.backfill_spliceai_max_ds(vav, chunk_size=chunk_size)
            logging.info("VariantAnnotationVersion %s: updated %d rows", vav.pk, updated)
            if not options["skip_index"]:
                self._create_partition_index(vav)
            vav.backfilled_spliceai_max_ds = True
            vav.save(update_fields=["backfilled_spliceai_max_ds"])

    @staticmethod
    def _create_partition_index(vav: VariantAnnotationVersion) -> None:
        base = VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION
        partition_table = vav.get_partition_table(base_table_name=base)
        index_name = f"{partition_table}_spliceai_max_ds_idx"
        sql = (
            f'CREATE INDEX CONCURRENTLY IF NOT EXISTS "{index_name}" '
            f'ON "{partition_table}" ("spliceai_max_ds")'
        )
        logging.info("VariantAnnotationVersion %s: %s", vav.pk, sql)
        with connection.cursor() as cursor:
            cursor.execute(sql)
