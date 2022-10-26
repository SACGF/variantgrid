from typing import List, Dict

from django.db.models import QuerySet

from classification.enums import EvidenceKeyValueType
from classification.models import EvidenceKey, Classification, ClassificationModification


class EvidenceKeyToUnit:

    def __init__(self, key_names: List[str]):
        self.e_keys: List[EvidenceKey] = []
        for key_name in key_names:
            e_key = EvidenceKey.objects.get(pk=key_name)
            self.e_keys.append(e_key)
            if e_key.value_type != EvidenceKeyValueType.UNIT:
                raise ValueError(f"EvidenceKey {e_key.key} is not a unit, it's {e_key.value_type}")

    def migrate_all(self):
        self.migrate(qs=Classification.objects.all())

    def migrate(self, qs: QuerySet, dry_run: bool = False):
        # changes: List[str] = []
        c: Classification
        for c in qs:
            converted_in_version: Dict[str, int] = {}
            cm: ClassificationModification
            for cm in ClassificationModification.objects.filter(classification__id=c.id).order_by('-id'):
                if evidence := cm.delta:
                    for e_key in self.e_keys:
                        if e_key_blob := evidence.get(e_key.key):
                            if note := e_key_blob.get('note'):
                                if 'Converted from' in note:
                                    converted_in_version[e_key.key] = cm.id

            for cm in ClassificationModification.objects.filter(classification__id=c.id).order_by('id'):
                raw_evidence = cm.published_evidence
                delta = cm.delta
                modified = False
                for e_key in self.e_keys:
                    if earliest_version := converted_in_version.get(e_key.key):
                        if earliest_version > cm.id:
                            for data_name, data in [('evidence', raw_evidence), ('delta', delta)]:
                                if data:
                                    if e_key_blob := data.get(e_key.key):
                                        if value := e_key_blob.get('value'):
                                            try:
                                                e_key_blob['value'] = float(value) / 100
                                                e_key_blob['note'] = f"Converted from '{value}'%"
                                                modified = True
                                                # changes.append(f"Version {c.id}.{cm.id}:{e_key.key} in {data_name} value to {float(value) / 100}")
                                            except:
                                                # changes.append(f"Version {c.id}.{cm.id}:{e_key.key} in {data_name} found a value that's not a float: {value}")
                                                print(f"Version {c.id}.{cm.id}:{e_key.key} in {data_name} found a value that's not a float: {value}")

                if not dry_run and modified:
                    cm.save(update_fields=['published_evidence', 'delta'])
                    print(f"Updated history of classification {c.id} version {cm.id}")
        # return changes
