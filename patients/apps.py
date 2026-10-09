from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from patients.signals import extraction_match_health_check  # noqa: F401
        from patients.signals.external_pk_search import search_external_pk
        from patients.signals.patient_search import patient_search
        from patients.signals.specimen_search import extraction_search, specimen_search
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(search_external_pk, patient_search, specimen_search, extraction_search)
