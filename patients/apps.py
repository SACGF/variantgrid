from django.apps import AppConfig


class PatientsConfig(AppConfig):
    name = 'patients'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        pass
        # pylint: enable=import-outside-toplevel,unused-import
