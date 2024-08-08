import glob
import os

from django.conf import settings
from django.core.management.base import BaseCommand

from analysis.analysis_import_export import analysis_import
from analysis.models import AnalysisTemplateType, AnalysisTemplate, AnalysisTemplateVersion
from library.guardian_utils import admin_bot, add_public_group_read_permission
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--replace', default=False, action='store_true')

    def handle(self, *args, **options):
        ANALYSIS_TEMPLATES_DIR = os.path.join(settings.BASE_DIR, "analysis", "data", "analysis_templates")

        user = admin_bot()
        genome_build = GenomeBuild.grch37()  # Doesn't matter for templates
        replace = options["replace"]

        for filename in glob.glob(f"{ANALYSIS_TEMPLATES_DIR}/*.json"):
            analysis = analysis_import(user, genome_build, filename, replace=replace)
            if analysis is None:
                print(f"Skipped {filename}...")
                continue
            else:
                print(f"Created: {filename}")

            analysis.template_type = AnalysisTemplateType.TEMPLATE
            analysis.save()
            add_public_group_read_permission(analysis)

            AnalysisTemplate.objects.filter(name=analysis.name).delete()  # Clear existing
            analysis_template = AnalysisTemplate.objects.create(name=analysis.name, user=user, analysis=analysis)

            analysis_snapshot = analysis.clone()
            analysis_snapshot.template_type = AnalysisTemplateType.SNAPSHOT
            analysis_snapshot.visible = False
            analysis_snapshot.save()
            add_public_group_read_permission(analysis_snapshot)

            analysis_name_template = "%(template)s for %(input)s"
            AnalysisTemplateVersion.objects.create(template=analysis_template,
                                                   version=1,
                                                   analysis_name_template=analysis_name_template,
                                                   analysis_snapshot=analysis_snapshot)

            print(f"Created template: {analysis_template}")
