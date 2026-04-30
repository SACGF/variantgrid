import logging

from django.core.management.base import BaseCommand
from django.db import connection

from annotation.models.models import VariantAnnotation, VariantAnnotationVersion


class Command(BaseCommand):
    help = "Backfill VariantAnnotation.spliceai_max_ds for every VariantAnnotationVersion partition. Idempotent."

    def add_arguments(self, parser):
        parser.add_argument("--vav-id", type=int,
                            help="Only process this VariantAnnotationVersion (default: all)")
        parser.add_argument("--skip-index", action="store_true",
                            help="Skip CREATE INDEX CONCURRENTLY step")

    def handle(self, *args, **options):
        qs = VariantAnnotationVersion.objects.all().order_by("pk")
        if vav_id := options.get("vav_id"):
            qs = qs.filter(pk=vav_id)

        for vav in qs:
            updated = VariantAnnotation.backfill_spliceai_max_ds(vav)
            logging.info("VariantAnnotationVersion %s: updated %d rows", vav.pk, updated)
            if not options["skip_index"]:
                self._create_partition_index(vav)

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
