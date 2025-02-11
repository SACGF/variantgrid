from django.core.management import BaseCommand

from annotation.regexes import DbRegexes, db_ref_regexes
from classification.enums import EvidenceKeyValueType
from classification.models import VCDataDict, EvidenceKeyMap, VCBlobKeys, Classification
from library.utils import JsonDiffs


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--update', action='store_true',
                            help="Updates the data")

    def handle(self, *args, **options):
        should_update = options["update"]
        _key_to_regex = {
            'db_rs_id': DbRegexes.SNP,
            'uniprot_id': DbRegexes.UNIPROTKB,
            'clinvar_variation_id': DbRegexes.CLINVAR,
            'cosmic_id': DbRegexes.COSMIC,
            'gtr_id': DbRegexes.GTR
        }

        def fix_evidence(evidence: VCDataDict) -> bool:
            had_changes = False
            for e_key in EvidenceKeyMap.cached().all_keys:
                if blob := evidence.get(e_key.key):
                    value = blob.get(VCBlobKeys.VALUE.value)
                    note = blob.get(VCBlobKeys.NOTE.value)
                    new_db_refs = []
                    if value and e_key.value_type in (EvidenceKeyValueType.FREE_ENTRY, EvidenceKeyValueType.TEXT_AREA):
                        scan_value = str(value)
                        if results := db_ref_regexes.search(scan_value, default_regex=_key_to_regex.get(e_key.key)):
                            new_db_refs += [r.to_json() for r in results]

                    if note:
                        if results := db_ref_regexes.search(note):
                            new_db_refs += [r.to_json() for r in results]
                    if not new_db_refs:
                        new_db_refs = None
                    old_db_refs = blob.get(VCBlobKeys.DB_REFS.value)
                    if diffs := JsonDiffs.differences(old_db_refs, new_db_refs):
                        had_changes = True
                        blob[VCBlobKeys.DB_REFS.value] = new_db_refs
            return had_changes

        change_count = 0
        for i, c in enumerate(Classification.objects.iterator()):
            has_changes = fix_evidence(c.evidence)
            if has_changes:
                if should_update:
                    c.save(update_fields=['evidence'])
                change_count += 1

            for cm in c.classificationmodification_set.filter(published=True).all():
                if fix_evidence(cm.published_evidence):
                    if should_update:
                        cm.save(update_fields=['published_evidence'])

            if i % 100 == 0:
                print(f"processing classification {i}")
        print(f"Changes = {change_count}")

