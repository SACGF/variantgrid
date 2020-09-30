import celery
from sync.models.models import SyncDestination


@celery.task
def sync_all():
    print('Running all syncs')
    for sync_dest in SyncDestination.objects.all():
        sync_dest.run()
