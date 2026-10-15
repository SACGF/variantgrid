from django.dispatch import receiver

from classification.models.classification_import_run import (
    ClassificationImportRun,
    classification_imports_complete_signal,
)
from classification.tasks.classification_citation_tasks import prefetch_classification_citations


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def prefetch_citations_after_import(**kwargs):
    """ Imported classifications export their citations by title only once prefetched, so do it now rather
        than waiting for the nightly run (#1687) """
    prefetch_classification_citations.delay()
