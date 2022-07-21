from django.apps import AppConfig


class OntologyConfig(AppConfig):
    name = 'eventlog'

    def ready(self):
        from eventlog.signals import active_users_health_check  # pylint: disable=unused-import
