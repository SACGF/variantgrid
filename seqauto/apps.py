

from django.apps import AppConfig


class SeqautoConfig(AppConfig):
    name = 'seqauto'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        pass
        # pylint: enable=import-outside-toplevel,unused-import
