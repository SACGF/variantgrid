from django.apps import AppConfig


# noinspection PyUnresolvedReferences
class ClassificationConfig(AppConfig):
    name = 'classification'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        import classification.signals  # noqa: F401  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel
