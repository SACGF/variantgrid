

from django.apps import AppConfig


class SeqautoConfig(AppConfig):
    name = 'seqauto'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from seqauto.signals import (  # noqa: F401
            enrichment_kit_search,
            experiment_search,
            seqauto_integration_status,
            sequencing_run_search,
        )
        # pylint: enable=import-outside-toplevel,unused-import
