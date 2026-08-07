from django.apps import AppConfig


class PedigreeConfig(AppConfig):
    name = 'pedigree'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        pass
        # pylint: enable=import-outside-toplevel,unused-import
