from cyvcf2 import Reader
from django.core.management.base import BaseCommand

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from library.genomics.vcf_utils import cyvcf2_header_types


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

        ALOFT = ["Aloft_Confidence", "Aloft_pred", "Aloft_prob_Dominant", "Aloft_prob_Recessive", "Aloft_prob_Tolerant"]

        count = 0
        num_found = 0
        for v in reader:
            count += 1
            csq = v.INFO.get("CSQ")
            found_in_variant = False
            found_in_transcript = False

            for transcript_csq in csq.split(","):
                td = dict(zip(vep_columns, transcript_csq.split("|")))
                for aloft_key in ALOFT:
                    if info_val := td[aloft_key]:
                        values = info_val.split("&")
                        if any([v and v != '.' for v in values]):
                            found_in_transcript = True
                            print(f"{aloft_key}={info_val=}")
                if found_in_transcript:
                    print("-----")
                    found_in_variant = True
                    found_in_transcript = False

            if found_in_variant:
                num_found += 1
                print("=" * 50)

                #if mane_select := td.get("MANE_SELECT"):
                #    print(f'mane_select={mane_select}, {td["HGVSc"]}')

                #if td["PICK"]:
                #    print(f'{td["SYMBOL"]} {td["HGVSc"]}, impact: {td["IMPACT"]}, revel: {td["REVEL_score"]}')

        print(f"{count=} variants, {num_found=}")
