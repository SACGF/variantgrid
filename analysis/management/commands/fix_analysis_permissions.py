import sys

from django.core.management.base import BaseCommand

from analysis.analysis_import_export import analysis_export_to_file
from analysis.models.models_analysis import Analysis

# Analysis templates didn't set permissions correctly until fix on October 19th 2020
# Need to fix legacy data created in systems before this fix: https://github.com/SACGF/variantgrid/issues/84
# Can be removed once environments fixed
from library.guardian_utils import assign_permission_to_user_and_groups


class Command(BaseCommand):
    def handle(self, *args, **options):
        for analysis in Analysis.objects.filter(analysistemplaterun__isnull=False):
            assign_permission_to_user_and_groups(analysis.user, analysis)
