import logging
import time

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.models import VariantAnnotationVersion, VariantAnnotation
from genes.hgvs import HGVSMatcher


class Command(BaseCommand):
    def handle(self, *args, **options):
        update_time = 5

        for vav in VariantAnnotationVersion.objects.all():
            matcher = HGVSMatcher(vav.genome_build)
            va_qs = VariantAnnotation.objects.filter(version=vav, hgvs_g__isnull=True)  # Only blank ones
            total_todo = va_qs.count()
            print(f"{vav.genome_build}/{vav.vep}: {total_todo} records to update...")
            if total_todo == 0:
                continue
            va_qs = va_qs.select_related("variant", "variant__locus", "variant__locus__contig",
                                         "variant__locus__ref", "variant__alt")
            records = []
            batch_size = 1000
            start = None
            processed = 0
            total_done = 0

            for va in va_qs.iterator():
                if start is None:
                    start = time.time()

                svlen = va.variant.svlen or 0
                if abs(svlen) <= settings.HGVS_MAX_SEQUENCE_LENGTH:
                    try:
                        va.hgvs_g = matcher.variant_to_g_hgvs(va.variant)
                        records.append(va)
                    except Exception as e:
                        logging.error("Skipped %s: %s", va.variant, e)
                        continue

                    if len(records) == batch_size:
                        self._update_variant_annotation(records)
                        processed += len(records)

                        now = time.time()
                        time_taken = now - start
                        if time_taken >= update_time:
                            total_done += processed
                            percent = 100 * total_done / total_todo

                            print(f"{percent:.1f}% complete... (last {time_taken:.1f} secs = {int(processed / time_taken)} records/sec)")
                            start = now
                            processed = 0

                        records = []

            if records:
                self._update_variant_annotation(records)

    @staticmethod
    def _update_variant_annotation(records: list[VariantAnnotation]):
        VariantAnnotation.objects.bulk_update(records, fields=["hgvs_g"])
