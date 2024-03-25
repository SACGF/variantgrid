import json
from collections import Counter, defaultdict

from cyvcf2 import Reader
from django.core.management.base import BaseCommand

from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from library.genomics.vcf_utils import cyvcf2_header_types


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', required=True)
        parser.add_argument('--columns', help="comma separated list of columns")
        parser.add_argument('--print-lines', action='store_true', help="Print info about each line")
        parser.add_argument('--pick', action='store_true', help="Print pick for each record")
        parser.add_argument('--json', action='store_true', help="Print records as JSON")

    def handle(self, *args, **options):
        vcf = options["vcf"]
        columns = options["columns"]
        print_lines = options["print_lines"]
        pick = options["pick"]
        do_json = options["json"]

        reader = Reader(vcf)

        header_types = cyvcf2_header_types(reader)
        infos = header_types["INFO"]
        vep_columns = BulkVEPVCFAnnotationInserter._get_vep_columns_from_csq(infos)
        if columns:
            for column in columns.split(","):
                if column not in vep_columns:
                    raise ValueError(f"{column=} not in {vep_columns=}")
        else:
            # Just use all columns
            columns = vep_columns

        field_values = defaultdict(Counter)

        count = 0
        num_found = 0
        for v in reader:
            count += 1
            csq = v.INFO.get("CSQ")
            found_in_variant = False
            found_in_transcript = False

            for transcript_csq in csq.split(","):
                td = dict(zip(vep_columns, transcript_csq.split("|")))
                if do_json:
                    j = json.dumps(td, indent=4)
                    print(j)
                else:
                    if pick:
                        if td.get("PICK"):
                            if variant_id := td.get("variant_id"):
                                print(f"{variant_id=}")
                            print(td)
                            print("-" * 50)

                    for nk in columns:
                        if info_val := td[nk]:
                            values = info_val.split("&")
                            values = [v for v in values if v and v != '.']
                            if any(values):
                                found_in_transcript = True
                                if print_lines:
                                    print(f"{nk}={info_val=}")
                                for v in values:
                                    field_values[nk][v] += 1

                    if found_in_transcript:
                        if print_lines:
                            print("-----")
                        found_in_variant = True
                        found_in_transcript = False

            if found_in_variant:
                num_found += 1
                if print_lines:
                    print("=" * 50)

                #if mane_select := td.get("MANE_SELECT"):
                #    print(f'mane_select={mane_select}, {td["HGVSc"]}')

                #if td["PICK"]:
                #    print(f'{td["SYMBOL"]} {td["HGVSc"]}, impact: {td["IMPACT"]}, revel: {td["REVEL_score"]}')

        if not do_json:
            print(f"{count=} variants, {num_found=}")

            for field, counts in field_values.items():
                print("---------------")
                print(f"{field=}")
                if len(counts) > 20:
                    print(f"counts = {len(counts)} distinct values")
                else:
                    for v, count in counts.most_common():
                        print(f"{v} = {count}")

            empty_columns = [column for column in columns if column not in field_values]
            if empty_columns:
                print("-" * 50)
                print(f"Empty columns: {','.join(empty_columns)}")
