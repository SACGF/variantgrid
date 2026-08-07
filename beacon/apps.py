from django.apps import AppConfig


class BeaconConfig(AppConfig):
    name = 'beacon'
    verbose_name = "Beacon v2"

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        pass
        # pylint: enable=import-outside-toplevel,unused-import
