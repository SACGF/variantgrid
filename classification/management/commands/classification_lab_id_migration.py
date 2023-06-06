from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import List, Dict

from django.core.management import BaseCommand
from classification.enums import SpecialEKeys
from classification.models import Classification
from genes.hgvs import CHGVS
from snpdb.models import Lab
import re


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
        base_id = self.genome_build_prefix + ":" + self.c_hgvs_obj.without_transcript_version.full_c_hgvs
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

        clash_count = 0
        clash_with_clin_sig_change_count = 0
        clash_with_two_plus_active = 0

        for new_id, old_ids in clash_dict.items():
            if len(old_ids) > 1:
                # old_id_list = ", ".join(old_ids)
                clin_sigs = set()
                for classification in old_ids:
                    clin_sigs.add(classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

                if clash_count == 0:
                    old_id_list = ", ".join(classification.lab_record_id for classification in old_ids)
                    print(f"e.g. new_id <- {old_id_list}")

                non_withdrawn = set()
                for classification in old_ids:
                    if not classification.withdrawn:
                        non_withdrawn.add(classification)
                if len(non_withdrawn) > 1:
                    clash_list = ", ".join(
                        f"{classification.lab_record_id} {classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)} withdrawn={classification.withdrawn}"
                        for classification in old_ids)
                    print(f"two+ non withdrawn {new_id} <- {clash_list}")
                    clash_with_two_plus_active += 1

                clash_count += 1
                if len(clin_sigs) > 1:
                    clash_list = ", ".join(f"{classification.lab_record_id} {classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)} withdrawn={classification.withdrawn}" for classification in old_ids)
                    print(f"clash {new_id} <- {clash_list}")
                    clash_with_clin_sig_change_count += 1

        print(f"Clash count = {clash_count}")
        print(f"Clash where 2+ aren't withdrawn = {clash_with_two_plus_active}")
        print(f"Clashes with clin sig change = {clash_with_clin_sig_change_count}")
