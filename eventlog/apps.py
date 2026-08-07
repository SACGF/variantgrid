from django.apps import AppConfig


class OntologyConfig(AppConfig):
    name = 'eventlog'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from eventlog.signals import (
            active_users_health_check,  # noqa: F401  # pylint: disable=unused-import
        )
        # pylint: enable=import-outside-toplevel
