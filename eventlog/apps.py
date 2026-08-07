from django.apps import AppConfig


class OntologyConfig(AppConfig):
    name = 'eventlog'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        pass  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel
