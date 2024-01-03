# DO NOT import any models, so keep this safe to be used by migrations
import contextlib
from typing import Optional


class BulkUpdator:

    def __init__(self, model, fields: list[str]):
        self.model = model
        self.fields = fields
        self.batch = []
        self.records_updated = 0

    def append(self, obj):
        self.batch.append(obj)
        if len(self.batch) >= 1000:
            self.model.objects.bulk_update(self.batch, fields=self.fields)
            self.records_updated += len(self.batch)
            self.batch = []
            print(f"Bulk updated {self.model} x {self.records_updated}")

    def finish(self):
        if self.batch:
            self.model.objects.bulk_update(self.batch, fields=self.fields)
            self.records_updated += len(self.batch)
            self.batch = []
        print(f"Bulk updated {self.model} x {self.records_updated}")

    @staticmethod
    @contextlib.contextmanager
    def instance(model, fields: list[str]):
        updator = BulkUpdator(model=model, fields=fields)
        try:
            yield updator
        finally:
            updator.finish()


class EvidenceKeyRenamer:
    """
    Allows you to rename an evidence key and it will rename all instances of it in variant classification evidence
    (current and historical in modifications) as well us migrating the evidence key
    """

    def __init__(self, apps, old_key: Optional[str] = None, new_key: Optional[str] = None, key_dict: Optional[dict[str, str]] = None):
        self.apps = apps
        self.EvidenceKey = apps.get_model("classification", "EvidenceKey")
        if not key_dict:
            key_dict = {}
        if old_key and new_key:
            key_dict[old_key] = new_key
        self.key_dict = key_dict
        if not key_dict:
            raise ValueError("Have nothing to rename")

    def _move_evidence_keys(self):
        for old_key, new_key in self.key_dict.items():
            if self.EvidenceKey.objects.filter(key=new_key).exists():
                raise ValueError(f"New evidence key {new_key} already exists")
            if old_e_key := self.EvidenceKey.objects.filter(key=old_key).first():
                print(f"Moving {old_key} to be {new_key}")
                old_e_key.key = new_key
                old_e_key.save()
                self.EvidenceKey.objects.filter(key=old_key).delete()
            else:
                raise ValueError(f"Old evidence key {old_key} doesn't exist")

    def _rename_in_dict(self, data: dict) -> bool:
        updated = False
        for old_key, new_key in self.key_dict.items():
            if data and old_key in data:
                data[new_key] = data[old_key]
                data.pop(old_key)
                updated = True
        return updated

    def run(self):
        Classification = self.apps.get_model("classification", "Classification")
        ClassificationModification = self.apps.get_model("classification", "ClassificationModification")

        self._move_evidence_keys()

        with BulkUpdator.instance(model=Classification, fields=["evidence"]) as batch_update:
            for vc in Classification.objects.iterator():
                if self._rename_in_dict(vc.evidence):
                    batch_update.append(vc)

        with BulkUpdator.instance(model=ClassificationModification, fields=["published_evidence", "delta"]) as batch_update:
            for vcm in ClassificationModification.objects.iterator():
                update_published = self._rename_in_dict(vcm.published_evidence)
                update_delta = self._rename_in_dict(vcm.delta)
                if update_published or update_delta:
                    batch_update.append(vcm)

    @staticmethod
    def rename(apps, old_key: Optional[str] = None, new_key: Optional[str] = None, key_dict: Optional[dict[str, str]] = None):
        EvidenceKeyRenamer(apps=apps, old_key=old_key, new_key=new_key, key_dict=key_dict).run()


class EvidenceSelectKeyRenamer:

    def __init__(self, apps, key: str, old_option: str, new_option: str):
        self.apps = apps
        self.EvidenceKey = apps.get_model("classification", "EvidenceKey")
        self.key = key
        self.old_option = old_option
        self.new_option = new_option

    def _update_evidence_key(self):
        ekey = self.EvidenceKey.objects.get(pk=self.key)
        options = ekey.options
        old_option_json = [option for option in options if option.get('key') == self.old_option]
        if old_option_json:
            old_option_json = old_option_json[0]
            old_option_json['key'] = self.new_option

        if self.old_option.lower() != self.new_option.lower():  # no need to alias if just different case
            option_json = [option for option in options if option.get('key') == self.new_option][0]

            aliases = option_json.get('aliases', [])
            if self.old_option not in aliases:
                aliases.append(self.old_option)
            option_json['aliases'] = aliases

        ekey.save()

    def _update_in_dict(self, data: dict):
        if data:
            blob = data.get(self.key)
            if not isinstance(blob, dict):
                return

            if value := blob.get('value'):
                if isinstance(value, list):
                    if self.old_option in value:
                        index = value.index(self.old_option)
                        value[index] = self.new_option
                elif value == self.old_option:
                    blob['value'] = self.new_option

    def _migrate_data(self):
        Classification = self.apps.get_model("classification", "Classification")
        ClassificationModification = self.apps.get_model("classification", "ClassificationModification")

        kwargs = {f"evidence__{self.key}__isnull": False}
        for vc in Classification.objects.filter(**kwargs):
            self._update_in_dict(vc.evidence)
            vc.save(update_fields=["evidence"])

        for vcm in ClassificationModification.objects.filter():
            self._update_in_dict(vcm.published_evidence)
            self._update_in_dict(vcm.delta)

            vcm.save(update_fields=["published_evidence", "delta"])

    def run(self):
        self._migrate_data()
        self._update_evidence_key()
