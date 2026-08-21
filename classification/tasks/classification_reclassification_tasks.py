import logging

import celery

from classification.models import ReclassificationEventBuilder


@celery.shared_task(queue='db_workers')
def reclassification_events_update():
    """
    Brings every clinical significance timeline up to date with the classifications published since the
    last run. Also the fallback the analytics page queues when a batch is too big to build in the request.
    @see ReclassificationEvent
    """
    result = ReclassificationEventBuilder.bring_up_to_date()
    logging.info("Reclassification timelines: rebuilt %d classifications (%d events), built to %s",
                 result.classifications, result.events, result.built_to)
