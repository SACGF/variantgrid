import logging

import celery

from classification.models import ReclassificationEventBuilder


@celery.shared_task(queue='db_workers')
def reclassification_events_reconcile():
    """
    Safety net for the classification_post_publish_signal receiver - rebuilds the timeline of any
    classification whose last published significance no longer matches the end of its timeline.
    @see ReclassificationEvent
    """
    needs_reconcile_qs = ReclassificationEventBuilder.classifications_needing_reconcile()
    classification_count = needs_reconcile_qs.count()
    if not classification_count:
        return

    rows_written = ReclassificationEventBuilder.rebuild(needs_reconcile_qs)
    logging.info("Reclassification reconcile rebuilt %d classifications (%d events)",
                 classification_count, rows_written)
