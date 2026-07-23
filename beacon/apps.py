from django.apps import AppConfig


class BeaconConfig(AppConfig):
    name = 'beacon'
    verbose_name = "Beacon v2"

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from beacon.signals import beacon_health_check
        # pylint: enable=import-outside-toplevel,unused-import
