

from django.apps import AppConfig


class SeqautoConfig(AppConfig):
    name = 'seqauto'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from seqauto.signals import experiment_search, sequencing_run_search
        # pylint: enable=import-outside-toplevel,unused-import
