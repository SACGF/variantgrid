#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from annotation.models import AnnotationRun, AnnotationRangeLock


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        if options["clear"]:
            logging.info("Clearing Everything")
            AnnotationRun.objects.all().update(annotation_range_lock=None)
            AnnotationRangeLock.objects.all().delete()
        else:
            logging.info("Clearing incomplete runs and locks")
            incomplete_locks_qs = AnnotationRangeLock.objects.filter(annotationrun__annotation_end__isnull=True)
            incomplete_annotation_runs_qs = AnnotationRun.objects.filter(annotation_range_lock__isnull=False, annotation_end__isnull=True)

            num_runs = incomplete_annotation_runs_qs.count()
            num_locks = incomplete_locks_qs.count()

            if num_runs or num_locks:
                incomplete_locks_qs.delete()
                incomplete_annotation_runs_qs.delete()
                logging.warning("Deleted %d incomplete annotation runs and %d locks", num_runs, num_locks)
            else:
                logging.warning("No incomplete runs or locks found")
