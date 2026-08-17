from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from patients.signals import (  # noqa: F401
            extraction_match_health_check,
            external_pk_search,
            patient_search,
            specimen_search,
        )
        # pylint: enable=import-outside-toplevel,unused-import
