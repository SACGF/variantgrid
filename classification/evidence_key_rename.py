# DO NOT import any models, so keep this safe to be used by migrations

class EvidenceKeyRenamer:
    """
    Allows you to rename an evidence key and it will rename all instances of it in variant classification evidence
    (current and historical in modifications) as well us migrating the evidence key
    """

    def __init__(self, apps, old_key: str, new_key: str):
        self.apps = apps
        self.EvidenceKey = apps.get_model("classification", "EvidenceKey")
        self.old_key = old_key
        self.new_key = new_key

    def _create_new_evidence_key(self):
        new_key = self.EvidenceKey.objects.filter(key=self.new_key).first()
        old_key = self.EvidenceKey.objects.filter(key=self.old_key).first()
        if not new_key and old_key:
            print(f"Copying {self.old_key} into new key {self.new_key}")
            new_key = old_key
            new_key.key = self.new_key
            new_key.save()
            return True
        print("Old key doesn't exist or new key already exists - not migrating")
        return False

    def _rename_in_dict(self, data: dict):
        if data and self.old_key in data:
            data[self.new_key] = data[self.old_key]
            data.pop(self.old_key)

    def _migrate_existing_evidence_key(self):
        Classification = self.apps.get_model("classification", "Classification")
        ClassificationModification = self.apps.get_model("classification", "ClassificationModification")

        kwargs = {f"evidence__{self.old_key}__isnull": False}
        for vc in Classification.objects.filter(**kwargs):
            self._rename_in_dict(vc.evidence)
            vc.save(update_fields=["evidence"])

        for vcm in ClassificationModification.objects.filter():
            self._rename_in_dict(vcm.published_evidence)
            self._rename_in_dict(vcm.delta)

            vcm.save(update_fields=["published_evidence", "delta"])

    def _delete_old_evidence_key(self):
        self.EvidenceKey.objects.filter(key=self.old_key).delete()

    def run(self):
        if self._create_new_evidence_key():
            self._migrate_existing_evidence_key()
            self._delete_old_evidence_key()

    @staticmethod
    def rename(apps, old_key: str, new_key: str):
        EvidenceKeyRenamer(apps=apps, old_key=old_key, new_key=new_key).run()


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
