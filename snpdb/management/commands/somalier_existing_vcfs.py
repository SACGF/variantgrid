from django.conf import settings
from django.core.management.base import BaseCommand

from snpdb.models import VCF, log_traceback, SomalierVCFExtract, GenomeBuild, SomalierRelatePairs
from snpdb.tasks.somalier_tasks import somalier_vcf_id, somalier_all_samples


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--genome-build")
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        if options.get("clear"):
            SomalierVCFExtract.objects.all().delete()
            SomalierRelatePairs.objects.all().delete()

        vcf_kwargs = {}
        if build_name := options.get("genome_build"):
            vcf_kwargs["genome_build"] = GenomeBuild.get_name_or_alias(build_name)

        if not settings.SOMALIER.get("enabled"):
            raise ValueError("settings.SOMALIER['enabled'] not enabled!")

        for vcf in VCF.objects.filter(somaliervcfextract__isnull=True, **vcf_kwargs):
            try:
                somalier_vcf_id(vcf.pk)
            except:
                log_traceback()

        somalier_all_samples()
