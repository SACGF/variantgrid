from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        from patients.signals import patient_search
        # pylint: enable=import-outside-toplevel
