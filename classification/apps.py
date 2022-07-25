from django.apps import AppConfig


class ClassificationConfig(AppConfig):
    name = 'classification'

    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        from classification.signals import classification_search  # pylint: disable=unused-import
        from classification.signals import classification_health_checks  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel
