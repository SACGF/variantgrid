

from django.apps import AppConfig


class SeqautoConfig(AppConfig):
    name = 'seqauto'

    # noinspection PyUnresolvedReferences
    def ready(self):
        from seqauto.signals import experiment_search
        from seqauto.signals import sequencing_run_search
        pass
