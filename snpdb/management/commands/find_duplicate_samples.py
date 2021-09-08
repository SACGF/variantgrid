#!/usr/bin/env python3

import logging
from collections import defaultdict

from django.core.management.base import BaseCommand

from snpdb.models import Sample


class Command(BaseCommand):
    def handle(self, *args, **options):
        duplicate_samples = defaultdict(list)
        qs = Sample.objects.filter(no_dna_control=False, samplestats__variant_count__isnull=False)
        for sample_id, name, vcf_name, variant_count in qs.values_list("pk", "name", "vcf__name", "samplestats__variant_count"):
            if "_NDC_" in name:
                continue
            key = f"{name}_{variant_count}"
            duplicate_samples[key].append((sample_id, name, vcf_name, variant_count))

        for dupes in duplicate_samples.values():
            if len(dupes) > 1:
                d = dupes[0]
                logging.info("%s (%d variants)", d[1], d[3])
                for d in dupes:
                    logging.info("%s (%s)", d[0], d[2])
                logging.info("*" * 10)
