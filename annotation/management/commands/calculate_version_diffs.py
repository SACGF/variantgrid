from django.core.management.base import BaseCommand
import logging

from annotation.models.models import EnsemblGeneAnnotationVersion
from annotation.models.models_version_diff import EnsemblGeneAnnotationVersionDiff, \
    VersionDiff


def calculate_version_diffs(version_klass, version_diff_klass):
    versions = list(version_klass.objects.all())
    version_from = None
    for version_to in versions:
        if version_from:
            try:
                version_diff = version_diff_klass.objects.get(version_from=version_from,
                                                              version_to=version_to)
            except version_diff_klass.DoesNotExist:
                version_diff = version_diff_klass(version_from=version_from,
                                                  version_to=version_to)
                version_diff.create_version_diffs()
        version_from = version_to


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        ANNOTATION_TYPES = [(EnsemblGeneAnnotationVersion, EnsemblGeneAnnotationVersionDiff)]

        if options.get("clear"):
            logging.info("Clearing Existing VersionDiffs")
            VersionDiff.objects.all().delete()

        for (version_klass, version_diff_klass) in ANNOTATION_TYPES:
            calculate_version_diffs(version_klass, version_diff_klass)
