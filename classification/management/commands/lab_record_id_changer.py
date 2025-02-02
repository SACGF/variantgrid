import re
from collections import defaultdict

from django.core.management import BaseCommand

from classification.enums import SpecialEKeys
from classification.models import Classification
from genes.hgvs import CHGVS
from snpdb.models import Lab, GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--lab_id', type=int, default=0)
        parser.add_argument('--commit', action='store_true')

    def handle(self, *args, **options):

        re_nucleotides = re.compile(r"(.*(?:del|dup))[ATCG]*$")

        lab = Lab.objects.get(pk=options['lab_id'])
        new_id_count = defaultdict(int)
        old_to_new = {}
        print(lab)
        pending_changes = []
        should_commit = options["commit"]
        for cr in Classification.objects.filter(lab=lab).iterator():
            c_hgvs_str = cr.get(SpecialEKeys.C_HGVS)
            c_hgvs = CHGVS(c_hgvs_str)
            genome_build_str = cr.get(SpecialEKeys.GENOME_BUILD)

            if c_hgvs.transcript and genome_build_str:
                genome_build = GenomeBuild.get_name_or_alias(genome_build_str)

                # new ID should exclude transcript version and gene symbol
                # only include genome_build if it's not GRCh38

                c_dot_normal = c_hgvs.c_dot
                if m := re_nucleotides.match(c_dot_normal):
                    c_dot_normal = m.group(1)
                c_dot_normal = c_dot_normal.replace("c.-", "c_minus")

                parts = [
                    genome_build.name,
                    c_hgvs.transcript_parts.identifier,
                    c_dot_normal
                ]

                record_id = "_".join(parts)
                record_id = record_id.replace('*', 'x')
                record_id = re.sub(r'\W+', '_', record_id)

                old_to_new[cr.lab_record_id] = record_id
                new_id_count[record_id] += 1

                if should_commit:
                    cr.lab_record_id = record_id
                    pending_changes.append(cr)
            else:
                print(f"Invalid genome build / c.HGVS {genome_build_str} / {c_hgvs_str}")

        for old_id, new_id in old_to_new.items():
            print(f"{old_id}\t{new_id}")

        is_valid = True
        for new_id, count in new_id_count.items():
            if count > 1:
                print(f"**** ID CLASH = {new_id}")
                is_valid = False

        if should_commit:
            if not is_valid:
                print("Can't update records due to ID Clash")
            else:
                Classification.objects.bulk_update(pending_changes, fields=["lab_record_id"])
