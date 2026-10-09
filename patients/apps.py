from importlib import import_module

from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    def import_models(self):
        super().import_models()
        # Not star-imported by patients.models: it imports ontology, which imports snpdb, which imports patients.models
        import_module(f"{self.name}.models.models_phenotype")

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from patients.signals import (  # noqa: F401
            ambiguous_acronym_denylist,
            external_pk_search,
            extraction_match_health_check,
            patient_search,
            phenotype_description,
            specimen_search,
        )
        # pylint: enable=import-outside-toplevel,unused-import
