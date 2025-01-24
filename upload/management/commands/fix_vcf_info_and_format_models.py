import logging
import os
from tempfile import NamedTemporaryFile

from cyvcf2 import Reader
from django.core.management.base import BaseCommand
from django.db.models.query_utils import Q

from library.genomics.vcf_utils import cyvcf2_header_types
from snpdb.models import VCF
from upload.vcf.vcf_import import create_vcf_info, create_vcf_format


class Command(BaseCommand):
    """
        See upload.vcf.bulk_genotype_vcf_processor.BulkGenotypeVCFProcessor
        v21 - started inserting info/format fields
        v25 - store VCFInfo/VCFFilter based on header

        So 21-24 we can just re-read the header and make these models. Prior to 21 we'll have to reload the VCF
        completely
    """
    def handle(self, *args, **options):
        importer_name = "PythonKnownVariantsImporter"
        v_info_format_min = 21
        v_info_format_max = 24

        q_header_has_info_or_format = Q(header__icontains='INFO') | Q(header__icontains='INFO')
        qs_old_version = VCF.objects.filter(q_header_has_info_or_format,
                                            uploadedvcf__vcf_importer__name=importer_name,
                                            uploadedvcf__vcf_importer__version__gte=v_info_format_min,
                                            uploadedvcf__vcf_importer__version__lte=v_info_format_max)

        logging.info("VCFs v%d - v%d with INFO/FORMAT header we can reprocess header only: %d",
                     v_info_format_min, v_info_format_max, qs_old_version.count())

        qs_can_add_models = qs_old_version.filter(vcfinfo__isnull=True, vcfformat__isnull=True)
        logging.info("Number without models: %d", qs_can_add_models.count())
        for vcf in qs_can_add_models:
            header_types = self.get_header_types(vcf)
            create_vcf_info(vcf, header_types.get("INFO", {}))
            create_vcf_format(vcf, header_types.get("FORMAT", {}))
            info_vals = vcf.vcfinfo_set.all().values_list("identifier", flat=True)
            format_vals = vcf.vcfformat_set.all().values_list("identifier", flat=True)
            logging.info("%s: INFO: %s, FORMAT: %s", vcf, ','.join(info_vals), ','.join(format_vals))


        qs_need_reload = VCF.objects.filter(q_header_has_info_or_format,
                                            uploadedvcf__vcf_importer__name=importer_name,
                                            uploadedvcf__vcf_importer__version__lt=v_info_format_min)
        if num_vcfs_need_reload := qs_need_reload.count():
            print(f"VCFs with INFO/FORMAT header that need to be reloaded {num_vcfs_need_reload}")

    @staticmethod
    def get_header_types(vcf):
        filename = vcf.uploadedvcf.uploaded_file.get_filename()
        if os.path.isfile(filename):
            reader = Reader(filename)
        else:
            # Make a temp file?
            with NamedTemporaryFile(mode='w+', delete=True) as temp_file:
                temp_file.write(vcf.header)
                temp_file.flush()  # Ensure the data is written to disk
                temp_file.seek(0)  # Move the pointer back to the start of the file
                reader = Reader(temp_file)

        header_types = cyvcf2_header_types(reader)
        return header_types