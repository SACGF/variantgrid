from django.apps import AppConfig

class ManualConfig(AppConfig):
    name = 'manual'

    def ready(self):
        pass
