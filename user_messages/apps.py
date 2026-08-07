from django.apps import AppConfig
from django.db.models.signals import post_save


class UserMessagesConfig(AppConfig):
    name = 'user_messages'
    verbose_name = "User Messages"

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from user_messages.models import Message
        from user_messages.signal_handlers import email_new_message_handler
        # pylint: enable=import-outside-toplevel,unused-import

        post_save.connect(email_new_message_handler, sender=Message)
