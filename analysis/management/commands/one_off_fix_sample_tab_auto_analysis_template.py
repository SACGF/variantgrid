"""
    This needs to be a management command so we can use code / model methods etc
"""

from django.conf import settings
from django.core.management import BaseCommand

from analysis.models import AnalysisTemplate
from analysis.models.nodes.node_utils import reload_analysis_nodes


class Command(BaseCommand):
    def handle(self, *args, **options):
        if analysis_template := AnalysisTemplate.objects.filter(name=settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE).first():
            analysis = analysis_template.analysis
            if sample_node := analysis.analysisnode_set.filter(
                    name='All variants').select_subclasses().first():
                if not all([sample_node.zygosity_ref, sample_node.zygosity_unk, sample_node.output_node]):
                    print(f"Updating {settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE} and creating new version...")
                    sample_node.zygosity_ref = True
                    sample_node.zygosity_unk = True
                    sample_node.output_node = True
                    sample_node.version += 1
                    sample_node.save()

                    reload_analysis_nodes(analysis.pk)

                    analysis_name_template = analysis_template.default_name_template()
                    analysis_template.new_version(analysis_name_template)

