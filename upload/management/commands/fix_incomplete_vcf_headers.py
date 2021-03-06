import os

from cyvcf2 import Reader
from django.core.management.base import BaseCommand
from django.db.models import Q

from snpdb.models import VCF
from upload.models import UploadedVCF


class Command(BaseCommand):
    def handle(self, *args, **options):
        q_no_header = Q(header__isnull=True)
        q_no_chrom_line = ~Q(header__icontains='#CHROM')

        num_fixed = 0
        missing_file = 0
        for vcf in VCF.objects.filter(q_no_header | q_no_chrom_line):
            try:
                filename = vcf.uploadedvcf.uploaded_file.get_filename()
                if os.path.exists(filename):
                    reader = Reader(filename)
                    vcf.header = reader.raw_header
                    vcf.save()
                    num_fixed += 1
                    if num_fixed % 100 == 0:
                        print(f"Fixed {num_fixed}")
                else:
                    missing_file += 1
            except UploadedVCF.DoesNotExist:
                missing_file += 1

        print(f"VCFs with incomplete headers - fixed: {num_fixed}, missing: {missing_file}")
