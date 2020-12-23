from django.conf import settings
from django.core.management.base import BaseCommand

from snpdb.models import VCF, log_traceback, SomalierVCFExtract
from snpdb.tasks.somalier_tasks import somalier_vcf_id, somalier_all_samples


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        if options.get("clear"):
            SomalierVCFExtract.objects.all().delete()

        if not settings.SOMALIER.get("enabled"):
            raise ValueError("settings.SOMALIER['enabled'] not enabled!")

        for vcf in VCF.objects.filter(somaliervcfextract__isnull=True):
            try:
                somalier_vcf_id(vcf.pk)
            except:
                log_traceback()

        somalier_all_samples()
