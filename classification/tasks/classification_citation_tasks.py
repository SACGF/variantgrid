"""
Nightly prefetch of the PubMed/PMC/Bookshelf citations referenced by published classifications (#1687), so
exports and the classification form read citations from the database instead of calling Entrez on demand.

Entry points: `prefetch_classification_citations` (celery beat) and `classification_citation_ids`.
"""
import celery
from django.db import connection

from annotation.models.models_citations import (
    CitationFetchRequest,
    CitationIdNormalized,
    CitationSource,
)
from library.utils import batch_iterator

ENTREZ_BATCH_SIZE = 100


def classification_citation_ids() -> set[CitationIdNormalized]:
    """ Every citation in the db_refs of a classification's last published evidence - what exports use """
    sql = """
        SELECT DISTINCT db_ref->>'db', db_ref->>'idx'
        FROM classification_classificationmodification cm,
            jsonb_path_query(cm.published_evidence, '$.*.db_refs[*]') db_ref
        WHERE cm.is_last_published
    """
    citation_ids = set()
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for db, idx in cursor.fetchall():
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
