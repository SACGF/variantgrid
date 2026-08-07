#!/usr/bin/env python3

import logging
from collections import defaultdict

from django.core.management.base import BaseCommand

from snpdb.models import CohortGenotypeStats


class Command(BaseCommand):
    def handle(self, *args, **options):
        duplicate_samples = defaultdict(list)
        # Per-sample CohortGenotypeStats rows (sample IS NOT NULL, filter_key NULL,
        # passing_filter=False) carry variant_count for each Sample.
        stats_qs = (CohortGenotypeStats.objects
                    .filter(sample__isnull=False, filter_key__isnull=True,
                            passing_filter=False, sample__no_dna_control=False,
                            variant_count__isnull=False)
                    .values_list("sample_id", "sample__name", "sample__vcf__name", "variant_count"))
        for sample_id, name, vcf_name, variant_count in stats_qs:
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
