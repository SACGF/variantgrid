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
            extraction_match_health_check,
            phenotype_description,
        )
        from patients.signals.external_pk_search import search_external_pk
        from patients.signals.patient_search import patient_search
        from patients.signals.specimen_search import extraction_search, specimen_search
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(search_external_pk, patient_search, specimen_search, extraction_search)
