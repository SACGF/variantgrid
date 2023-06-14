import re
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import List, Dict

from django.core.management import BaseCommand

from classification.enums import SpecialEKeys
from classification.models import Classification
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from snpdb.models import Lab


@dataclass
class LabRecordIdCalculator:
    genome_build: str
    c_hgvs: str
    lab_record_id: str

    @staticmethod
    def from_classification(classification: Classification) -> 'LabRecordIdCalculator':
        genome_build = classification.get(SpecialEKeys.GENOME_BUILD)
        c_hgvs = classification.get(SpecialEKeys.C_HGVS)
        lab_record_id = classification.lab_record_id
        return LabRecordIdCalculator(
            genome_build=genome_build,
            c_hgvs=c_hgvs,
            lab_record_id=lab_record_id
        )

    @cached_property
    def c_hgvs_obj(self):
        return CHGVS(self.c_hgvs)

    @cached_property
    def genome_build_prefix(self):
        if '37' in self.genome_build or '19' in self.genome_build:
            return "hg19"
        else:
            return "hg38"

    @cached_property
    def derive_new_id(self):
        base_id = self.genome_build_prefix + ":" + self.c_hgvs_obj.without_transcript_version.without_gene_symbol_str
        base_id = base_id.replace("*", "x")
        base_id = re.sub(r'\W+', '_', base_id)
        return base_id


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--lab', type=str, required=True)

    def handle(self, *args, **options):
        lab_id = options["lab"]
        lab = Lab.objects.get(group_name=lab_id)
        clash_dict: Dict[str, List[Classification]] = defaultdict(list)

        classification: Classification
        for classification in Classification.objects.filter(lab=lab).iterator():
            lab_id_calc = LabRecordIdCalculator.from_classification(classification)
            new_lab_record_id = lab_id_calc.derive_new_id
            clash_dict[new_lab_record_id].append(classification)

        # clash_count = 0
        # clash_with_clin_sig_change_count = 0
        # clash_with_two_plus_active = 0

        print("Operation\tLab Record ID\tc.HGVS\tTranscript Version\tGene Symbol\tClin Sig\tLast Imported")
        for new_id, records in clash_dict.items():
            records = [c for c in records if not c.withdrawn]

            if len(records) > 1:
                # old_id_list = ", ".join(records)

                def detailed_classification_str(c: Classification):
                    last_imported = c.created
                    if import_run := c.last_import_run:
                        last_imported = import_run.created
                    c_hgvs = CHGVS(c.imported_c_hgvs)
                    return f"{c.lab_record_id}\t{c_hgvs.full_c_hgvs}\t{c_hgvs.transcript_parts.version}\t{c_hgvs.gene_symbol}\t{c.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)}\t{last_imported:%Y-%m-%d}"

                records = sorted(records, key=lambda c: (CHGVS(c.imported_c_hgvs).transcript_parts.version, c.last_import_run.created), reverse=True)

                best = records[0]
                print(f"Keep\t{detailed_classification_str(best)}")

                for old in records[1:]:
                    if old.last_import_run.created > records[0].last_import_run.created:
                        raise ValueError("Smaller transcript was imported later than newer transcript")
                    print(f"Withdraw\t{detailed_classification_str(old)}")
                    # old.set_withdrawn(user=admin_bot(), withdraw=True)
                    # print(f"Withdrawing {old.lab_record_id} in favour of {records[0].lab_record_id}")
        #
        #         clin_sigs = set()
        #         for classification in records:
        #             clin_sigs.add(classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
        #
        #         if clash_count == 0:
        #             old_id_list = ", ".join(classification.lab_record_id for classification in records)
        #             print(f"e.g. new_id <- {old_id_list}")
        #
        #         non_withdrawn = set()
        #         for classification in records:
        #             if not classification.withdrawn:
        #                 non_withdrawn.add(classification)
        #         if len(non_withdrawn) > 1:
        #             clash_list = ", ".join(
        #                 f"{classification.lab_record_id} {classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)} withdrawn={classification.withdrawn}"
        #                 for classification in records)
        #             print(f"two+ non withdrawn {new_id} <- {clash_list}")
        #             clash_with_two_plus_active += 1
        #
        #         clash_count += 1
        #         if len(clin_sigs) > 1:
        #             clash_list = ", ".join(f"{classification.lab_record_id} {classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)} withdrawn={classification.withdrawn}" for classification in records)
        #             print(f"clash {new_id} <- {clash_list}")
        #             clash_with_clin_sig_change_count += 1
        #
        # print(f"Clash count = {clash_count}")
        # print(f"Clash where 2+ aren't withdrawn = {clash_with_two_plus_active}")
        # print(f"Clashes with clin sig change = {clash_with_clin_sig_change_count}")
