from django.apps import AppConfig


class SyncConfig(AppConfig):
    name = 'sync'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        pass
        # pylint: enable=import-outside-toplevel,unused-import
