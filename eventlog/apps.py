from django.apps import AppConfig


class OntologyConfig(AppConfig):
    name = 'eventlog'

    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        from eventlog.signals import active_users_health_check  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel
