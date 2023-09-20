from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from patients.signals import external_pk_search, patient_search
        # pylint: enable=import-outside-toplevel,unused-import
