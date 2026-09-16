from django.core.management.base import BaseCommand

from library.guardian_utils import admin_bot
from snpdb.models import AlleleConversionTool, GenomeBuild
from snpdb.tasks.liftover_tasks import liftover_alleles


class Command(BaseCommand):
    category = "ops"
    help = "Lifts over any alleles not in both genome builds"

    def add_arguments(self, parser):
        parser.add_argument('--retry-tool', choices=[act.value for act in AlleleConversionTool],
                            help="Only re-attempt this conversion tool, for the alleles it has already failed on")

    def handle(self, *args, **options):
        user = admin_bot()
        retry_tool = options["retry_tool"]
        for genome_build in GenomeBuild.builds_with_annotation():
            liftover_alleles(user.username, genome_build.name, retry_tool)
