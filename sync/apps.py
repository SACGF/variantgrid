from django.apps import AppConfig


class SyncConfig(AppConfig):
    name = 'sync'

    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from sync.signals import sync_health_check
        # pylint: enable=import-outside-toplevel,unused-import
