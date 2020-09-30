#!/usr/bin/env python3

from django.core.management.base import BaseCommand

from annotation.models import SampleVariantAnnotationStats, \
    SampleVariantAnnotationStatsPassingFilter, SampleEnsemblGeneAnnotationStats, \
    SampleEnsemblGeneAnnotationStatsPassingFilter, SampleClinVarAnnotationStats, \
    SampleClinVarAnnotationStatsPassingFilter
from annotation.tasks.calculate_sample_stats import calculate_needed_stats
from snpdb.models import SampleStats, SampleStatsPassingFilter


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--async', action='store_true')
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        if options.get("clear"):
            print("Clearing all sample stats")
            CLASSES = [SampleStats,
                       SampleStatsPassingFilter,
                       SampleVariantAnnotationStats,
                       SampleVariantAnnotationStatsPassingFilter,
                       SampleEnsemblGeneAnnotationStats,
                       SampleEnsemblGeneAnnotationStatsPassingFilter,
                       SampleClinVarAnnotationStats,
                       SampleClinVarAnnotationStatsPassingFilter]

            for clazz in CLASSES:
                clazz.objects.all().delete()

        run_async = options['async']

        calculate_needed_stats(run_async)
