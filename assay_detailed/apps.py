from django.apps import AppConfig


# noinspection PyUnresolvedReferences
class AssayDetailedConfig(AppConfig):
    name = 'assay_detailed'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        pass
        # pylint: enable=import-outside-toplevel
