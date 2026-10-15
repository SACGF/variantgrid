from django.core.cache import cache
from django.dispatch import receiver

from classification.models.classification_import_run import (
    ClassificationImportRun,
    classification_imports_complete_signal,
)
from classification.tasks.classification_citation_tasks import prefetch_classification_citations

PREFETCH_DELAY_SECONDS = 5 * 60
_PREFETCH_QUEUED_KEY = "prefetch-classification-citations-queued"


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def prefetch_citations_after_import(**kwargs):
    """ Imported classifications export their citations by title only once prefetched, so do it soon rather than
        waiting for the nightly run (#1687). Imports arrive in bursts (one run per lab in a sync), so the first one
        queues a single delayed run that scans everything published by the time it starts; concurrent runs would
        fetch the same new citations from Entrez twice. The nightly run catches anything a burst straddles. """
    if cache.add(_PREFETCH_QUEUED_KEY, True, PREFETCH_DELAY_SECONDS):
        prefetch_classification_citations.apply_async(countdown=PREFETCH_DELAY_SECONDS)
