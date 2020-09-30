from django.core.management.base import BaseCommand

from snpdb.models.models_vcf import VCF
from snpdb.models.vcf_replace_data import vcf_replace_data


class Command(BaseCommand):
    help = 'Copy data from a one VCF to another (with the same samples)'

    def add_arguments(self, parser):
        parser.add_argument('--old-vcf-id', type=int, required=True)
        parser.add_argument('--new-vcf-id', type=int, required=True)

    def handle(self, *args, **options):
        old_vcf_id = options.get("old_vcf_id")
        new_vcf_id = options.get("new_vcf_id")

        old_vcf = VCF.objects.get(pk=old_vcf_id)
        new_vcf = VCF.objects.get(pk=new_vcf_id)

        vcf_replace_data(old_vcf, new_vcf)
