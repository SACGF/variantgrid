import celery

from sync.models.models import SyncDestination


@celery.task
def sync_all():
    for sync_dest in SyncDestination.objects.filter(enabled=True):
        sync_dest.run()
