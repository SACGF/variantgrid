from django.apps import AppConfig


class ClassificationConfig(AppConfig):
    name = 'classification'

    def ready(self):
        from classification.signals import classification_search
        from classification.signals import classification_health_checks