import logging

import celery

from sync.models.models import SyncDestination


@celery.shared_task
def sync_all():
    destinations = list(SyncDestination.objects.filter(enabled=True))
    if not destinations:
        logging.info("sync_all: no enabled SyncDestinations")
        return

    logging.info("sync_all: %d enabled destination(s): %s",
                 len(destinations), ", ".join(str(sd) for sd in destinations))
    for sync_dest in destinations:
        sync_dest.run()
