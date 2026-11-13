from library.utils import first


class EvidenceKeyOptionUpdator:
    """
    Used in migration scripts to link a value to export to ClinVar to an evidence key option
    """

    def __init__(self, apps, key: str):
        self.apps = apps
        EvidenceKey = apps.get_model("classification", "EvidenceKey")
        self.evidence_key = EvidenceKey.objects.get(key=key)
        self.options: list[dict] = self.evidence_key.options
        if not isinstance(self.options, list):
            raise ValueError(f"EvidenceKey {key} does not have any options")

    def set_option(self, option_key: str, option_data_key: str, option_data_value: str):
        if option := first(option for option in self.options if option.get('key') == option_key):
            option[option_data_key] = option_data_value
        else:
            raise ValueError(f"No such option {self.evidence_key.key}.{option}")
        self.evidence_key.save()

    def upsert_option(self, option_data: dict):
        if option_key := option_data.get("key"):
            if option := first(option for option in self.options if option.get('key') == option_key):
                for key, value in option_data.items():
                    option[key] = value
            else:
                self.options.append(option_data)
        else:
            raise ValueError(f"Option must have a value for key")
        self.evidence_key.save()

    def set_clinvar_option(self, option_key: str, clinvar: str):
        self.set_option(option_key, "clinvar", clinvar)

    def set_testing_context_bucket(self, option_key: str, testing_context_bucket: str):
        self.set_option(option_key, "testing_context_bucket", testing_context_bucket)
