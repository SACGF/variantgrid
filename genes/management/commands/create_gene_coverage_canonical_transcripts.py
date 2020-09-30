#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from genes.models import GeneCoverageCanonicalTranscript, GeneCoverageCollection
from genes.tasks.gene_coverage_tasks import create_canonical_gene_coverage_for_enrichment_kit
from library.django_utils import get_field_counts
from seqauto.models import EnrichmentKit


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true', required=False)

    def handle(self, *args, **options):
        if options.get("clear"):
            GeneCoverageCanonicalTranscript.objects.all().delete()
        elif GeneCoverageCanonicalTranscript.objects.exists():
            logging.warning("Warning, there are already GeneCoverageCanonicalTranscript objects. You probably don't want to run this... Or maybe you want to run --clear")

        qs = GeneCoverageCollection.objects.exclude(genecoveragecanonicaltranscript__isnull=False)
        field_counts = get_field_counts(qs, "qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit__name")
        logging.info("Work to do:")
        logging.info(field_counts)

        for enrichment_kit in EnrichmentKit.objects.all():
            task = create_canonical_gene_coverage_for_enrichment_kit.si(enrichment_kit.pk)  # @UndefinedVariable
            task.apply_async()

        logging.info("Creating canonical coverage for GeneCoverage without enrichment kits (uploaded data)")
        task = create_canonical_gene_coverage_for_enrichment_kit.si(0)  # @UndefinedVariable
        task.apply_async()
