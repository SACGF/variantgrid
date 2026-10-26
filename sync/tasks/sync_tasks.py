import logging
from datetime import timedelta

import celery
from django.utils import timezone

from library.log_utils import report_exc_info
from sync.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.sync_run import run_sync

# A run still IN_PROGRESS younger than this is taken to be live, so the hourly sync_all leaves its destination alone.
# Older ones are from a worker that died without reaching run_sync's finally, and don't block.
IN_PROGRESS_RUN_MAX_AGE = timedelta(hours=12)


def _has_live_run(sync_destination: SyncDestination) -> bool:
    return SyncRun.objects.filter(destination=sync_destination,
                                  status=SyncStatus.IN_PROGRESS,
                                  created__gte=timezone.now() - IN_PROGRESS_RUN_MAX_AGE).exists()


@celery.shared_task
def sync_all():
    destinations = list(SyncDestination.objects.filter(enabled=True))
    if not destinations:
        logging.info("sync_all: no enabled SyncDestinations")
        return

    logging.info("sync_all: %d enabled destination(s): %s",
                 len(destinations), ", ".join(str(sd) for sd in destinations))
    failed = []
    for sync_dest in destinations:
        if _has_live_run(sync_dest):
            logging.info("sync_all: '%s' still has a run in progress, skipping", sync_dest)
            continue
        try:
            run_sync(sync_dest)
        except Exception:
            # run_sync has already marked the run FAILED; carry on with the other destinations
            report_exc_info({"sync_destination": str(sync_dest)})
            failed.append(str(sync_dest))

    if failed:
        raise RuntimeError(f"sync_all: failed for {', '.join(failed)}")
