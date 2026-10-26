

from django.apps import AppConfig

from library.django_utils.change_tracking import track_deletes


class SeqautoConfig(AppConfig):
    name = 'seqauto'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from seqauto.models.models_upload import (
            BackendVCF,
            UploadedDragenTSO500CombinedVariantOutput,
            UploadedDragenTSO500MetricsOutput,
        )
        from seqauto.signals import seqauto_integration_status, upload_signal_handlers  # noqa: F401
        from seqauto.signals.enrichment_kit_search import enrichment_kit_search
        from seqauto.signals.experiment_search import experiment_search
        from seqauto.signals.sequencing_run_search import sequencing_run_search
        from seqauto.signals.tso500_pair_search import tso500_pair_search
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel,unused-import

        search_registry.register(enrichment_kit_search, experiment_search, sequencing_run_search, tso500_pair_search)
        # The API Log counts deletes from the seqauto API (#1413) - the upload side's rows go with upload's
        upload_models = {BackendVCF, UploadedDragenTSO500CombinedVariantOutput, UploadedDragenTSO500MetricsOutput}
        track_deletes(*(model for model in self.get_models() if model not in upload_models))
