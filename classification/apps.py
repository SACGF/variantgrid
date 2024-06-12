from django.apps import AppConfig


# noinspection PyUnresolvedReferences
class ClassificationConfig(AppConfig):
    name = 'classification'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        import classification.signals  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel
