from django.apps import AppConfig


class BeaconConfig(AppConfig):
    name = 'beacon'
    verbose_name = "Beacon v2"

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from beacon.signals import beacon_health_check  # noqa: F401
        # pylint: enable=import-outside-toplevel,unused-import
