#!/usr/bin/env python3

import logging
from collections import Counter

from django.core.management.base import BaseCommand
from django.db.models import Q, Func, Value, F

from annotation.models import VariantAnnotation, VariantTranscriptAnnotation, VariantAnnotationVersion, \
    TranscriptVersion, defaultdict, AnnotationRun, VariantAnnotationPipelineType, AnnotationStatus
from genes.hgvs import HGVSMatcher
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    """ Only needs to be run on legacy systems that imported SV annotations before 2025-01-14 """

    def handle(self, *args, **options):
        for genome_build in GenomeBuild.builds_with_annotation():
            # We should have all the transcripts locally, don't fall back on ClinGen if we error as it'll
            # probably just error there too
            hgvs_matcher = HGVSMatcher(genome_build, clingen_resolution=False, allow_alternative_transcript_version=False)

            run_qs = AnnotationRun.objects.filter(pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                                                  status=AnnotationStatus.FINISHED)
            total = run_qs.count()
            # Pull down the
            for i, ar in enumerate(run_qs):
                for klass in [VariantAnnotation, VariantTranscriptAnnotation]:
                    hgvs_c_results = Counter()
                    qs = klass.objects.filter(annotation_run=ar, transcript_version__isnull=False)
                    qs = qs.exclude(hgvs_c=VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE)

                    variant_ids = []
                    transcript_accessions = []
                    records = []

                    for variant_id, pk, transcript_id, transcript_version in qs.values_list("variant_id", "pk", "transcript_version__transcript_id", "transcript_version__version"):
                        variant_ids.append(variant_id)
                        records.append(klass(pk=pk))
                        transcript_accession = TranscriptVersion.get_accession(transcript_id, transcript_version)
                        transcript_accessions.append(transcript_accession)

                    variant_coordinates_by_id = {}
                    for variant in Variant.objects.filter(pk__in=variant_ids).select_related("locus", "locus__ref", "alt"):
                        variant_coordinates_by_id[variant.pk] = variant.coordinate

                    for variant_id, transcript_accession, record in zip(variant_ids, transcript_accessions, records):
                        variant_coordinate = variant_coordinates_by_id[variant_id]
                        try:
                            hgvs_c = hgvs_matcher.variant_coordinate_to_hgvs_variant(variant_coordinate,
                                                                                     transcript_accession)
                            hgvs_c_results["ok"] += 1
                        except Exception as e:
                            hgvs_c = VariantAnnotation.SV_HGVS_ERROR_MESSAGE
                            hgvs_c_results["error"] += 1

                        record.hgvs_c = hgvs_c

                    if records:
                        logging.info("hgvs_c_results: %s", hgvs_c_results)
                        logging.info("%s - Bulk updating %d records", ar, len(records))
                        klass.objects.bulk_update(records, ["hgvs_c"], batch_size=1000)

                if i and (i % 10 == 0):
                    logging.info("Processed %d of %d records (%f)", i, total, 100*i/total)
