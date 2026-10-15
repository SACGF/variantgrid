"""
Prefetch of the PubMed/PMC/Bookshelf citations referenced by published classifications (#1687), so exports read
citations from the database instead of calling Entrez mid-download. Runs nightly (celery beat) for the refresh of
stale or errored citations, and after each classification import so newly imported ones are loaded straight away.

Entry points: `prefetch_classification_citations` (celery task) and `classification_citation_ids`.
"""
import celery

from annotation.models.models_citations import (
    CitationFetchRequest,
    CitationIdNormalized,
    CitationSource,
)
from classification.models.classification import ClassificationModification
from library.utils import batch_iterator

ENTREZ_BATCH_SIZE = 100


def classification_citation_ids() -> set[CitationIdNormalized]:
    """ Every citation in the db_refs of a classification's last published evidence - what exports use """
    citation_ids = set()
    modifications = ClassificationModification.objects.filter(is_last_published=True).only('published_evidence')
    for modification in modifications.iterator(chunk_size=2000):
        for db_ref in modification.db_refs:
            db, idx = db_ref.get('db'), db_ref.get('idx')
            if db and idx and CitationSource.from_legacy_code(db):
                try:
                    citation_ids.add(CitationIdNormalized.from_parts(source=db, index=idx))
                except ValueError:  # unparseable index typed in by a user
                    pass
    return citation_ids


@celery.shared_task
def prefetch_classification_citations():
    """ Loads citations never fetched, due a refresh or in error. Already loaded ones cost a DB query per batch """
    for batch in batch_iterator(sorted(classification_citation_ids()), batch_size=ENTREZ_BATCH_SIZE):
        CitationFetchRequest.fetch_all_now(batch)
