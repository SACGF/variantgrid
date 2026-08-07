from django.apps import AppConfig
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured


class MMEConfig(AppConfig):
    name = 'mme'
    verbose_name = "MatchMaker Exchange"

    def ready(self):
        """ EmailLog.send_mail is quietly a no-op when from_email is falsy or SEND_EMAILS
            is False. MME obliges us to notify depositors of matches, so refuse to start
            half-configured rather than find out from an email that never went. """
        if settings.MME_ENABLED and not (settings.MME_FROM_EMAIL and settings.SEND_EMAILS):
            raise ImproperlyConfigured(
                "MME_ENABLED requires MME_FROM_EMAIL and SEND_EMAILS to be set - "
                "MatchMaker Exchange must be able to notify depositors of matches")
