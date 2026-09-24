from django.apps import AppConfig


class SyncConfig(AppConfig):
    name = 'sync'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers and @register_sync_runner classes on import - noqa: F401 keeps the
        # unused-import autofix from silently unregistering them
        import sync.shariant  # noqa: F401
        from sync.signals import (
            sync_health_check,  # noqa: F401
            sync_integration_status,  # noqa: F401
        )
        # pylint: enable=import-outside-toplevel,unused-import
