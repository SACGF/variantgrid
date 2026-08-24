from django.apps import AppConfig


class VariantopediaConfig(AppConfig):
    name = 'variantopedia'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from variantopedia.signals import integration_health_check  # noqa: F401
        # pylint: enable=import-outside-toplevel,unused-import
