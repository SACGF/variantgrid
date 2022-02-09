from cyvcf2 import Reader
from django.core.management.base import BaseCommand

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from library.vcf_utils import cyvcf2_header_types


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', required=True)

    def handle(self, *args, **options):
        vcf = options["vcf"]
        reader = Reader(vcf)

        header_types = cyvcf2_header_types(reader)
        infos = header_types["INFO"]
        vep_columns = BulkVEPVCFAnnotationInserter._get_vep_columns_from_csq(infos)
        #print(vep_columns)

        for v in reader:
            #print(v)
            csq = v.INFO.get("CSQ")

            for transcript_csq in csq.split(","):
                td = dict(zip(vep_columns, transcript_csq.split("|")))
                # print(td)
                if mane_select := td.get("MANE_SELECT"):
                    print(f'mane_select={mane_select}, {td["HGVSc"]}')

                #if td["PICK"]:
                #    print(f'{td["SYMBOL"]} {td["HGVSc"]}, impact: {td["IMPACT"]}, revel: {td["REVEL_score"]}')