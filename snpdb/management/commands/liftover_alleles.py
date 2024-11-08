from django.core.management.base import BaseCommand

from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild
from snpdb.tasks.liftover_tasks import liftover_alleles


class Command(BaseCommand):
    help = "Lifts over any alleles not in both genome builds"

    def handle(self, **options):
        user = admin_bot()
        for genome_build in GenomeBuild.builds_with_annotation():
            liftover_alleles(user.username, genome_build.name)
