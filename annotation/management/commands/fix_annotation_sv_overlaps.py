#!/usr/bin/env python3

import logging
from collections import Counter

from django.core.management.base import BaseCommand
from django.db import transaction

from annotation.models import VariantAnnotation, VariantAnnotationVersion, VariantGeneOverlap, VEPSkippedReason
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import SVGeneOverlapResolver


def fix_annotation_sv_overlaps():
    """ Backfill overlapping_symbols + VariantGeneOverlap rows for long SVs (TOO_LONG)
        that VEP skipped. Idempotent: skips rows that already have overlapping_symbols set. """
    for vav in VariantAnnotationVersion.objects.filter(gene_annotation_release__isnull=False):
        logging.info("Processing %s", vav)
        resolver = SVGeneOverlapResolver(vav)

        qs = VariantAnnotation.objects.filter(
            version=vav,
            vep_skipped_reason=VEPSkippedReason.TOO_LONG,
            overlapping_symbols__isnull=True,
        ).select_related("variant__locus__contig", "variant__locus__ref", "variant__alt")

        results = Counter()
        total = qs.count()
        if not total:
            logging.info("Nothing to do for %s", vav)
            continue

        logging.info("%s - %d candidate variants", vav, total)
        for i, va in enumerate(qs.iterator()):
            variant = va.variant
            variant_coordinate = variant.coordinate
            symbols, gene_ids = resolver.get_overlaps(variant_coordinate)

            with transaction.atomic():
                if symbols:
                    va.overlapping_symbols = ",".join(sorted(symbols))
                    va.save(update_fields=["overlapping_symbols"])
                    results["with_overlaps"] += 1
                else:
                    results["no_overlaps"] += 1

                overlap_records = [
                    VariantGeneOverlap(
                        version=vav,
                        annotation_run=va.annotation_run,
                        variant=variant,
                        gene_id=gene_id,
                    )
                    for gene_id in gene_ids
                ]
                if overlap_records:
                    VariantGeneOverlap.objects.bulk_create(overlap_records, ignore_conflicts=True)
                    results["gene_overlap_rows"] += len(overlap_records)

            if i and (i % 100 == 0):
                logging.info("Processed %d of %d (%.1f%%)", i, total, 100 * i / total)

        logging.info("%s done: %s", vav, dict(results))


class Command(BaseCommand):
    """ Backfill long-SV gene overlaps (issue #1271).

        Long SVs (abs(svlen) > settings.ANNOTATION_VEP_SV_MAX_SIZE) are skipped by VEP and end up
        with vep_skipped_reason=TOO_LONG and no overlapping_symbols / VariantGeneOverlap rows.
        This command resolves gene overlaps locally from each VariantAnnotationVersion's
        gene_annotation_release transcripts and fills them in. """

    def handle(self, *args, **options):
        fix_annotation_sv_overlaps()
