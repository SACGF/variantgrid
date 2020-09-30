import sys

from django.core.management.base import BaseCommand

from analysis.analysis_import_export import analysis_export_to_file
from analysis.models.models_analysis import Analysis


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--analysis-id', type=int, required=True, help='Analysis PK')

    def handle(self, *args, **options):
        analysis_id = options["analysis_id"]

        analysis = Analysis.objects.get(pk=analysis_id)
        analysis_export_to_file(analysis, sys.stdout)