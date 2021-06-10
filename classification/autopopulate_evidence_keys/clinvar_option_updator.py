class ClinVarOptionUpdator:

    def __init__(self, apps, key: str):
        self.apps = apps
        EvidenceKey = apps.get_model("classification", "EvidenceKey")
        self.evidence_key = EvidenceKey.objects.get(key=key)
        self.options = self.evidence_key.options
        if not isinstance(self.options, list):
            raise ValueError(f"EvidenceKey {key} does not have any options")

    def set_clinvar_option(self, option_key: str, clinvar: str):
        if option := next(option for option in self.options if option.get('key') == option_key):
            option["clinvar"] = clinvar
        else:
            raise ValueError(f"No such option {self.key}.{option}")
        self.evidence_key.save()
