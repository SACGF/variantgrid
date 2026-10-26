from unittest.mock import patch

from django.test import TestCase

from sync.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.tasks.sync_tasks import sync_all


class SyncAllTestCase(TestCase):

    def setUp(self):
        config = {"type": "shariant", "direction": "upload", "sync_details": "test_shariant"}
        self.failing = SyncDestination.objects.create(name="Failing", config=config, enabled=True)
        self.working = SyncDestination.objects.create(name="Working", config=config, enabled=True)
        self.busy = SyncDestination.objects.create(name="Busy", config=config, enabled=True)
        SyncRun.objects.create(destination=self.busy, status=SyncStatus.IN_PROGRESS)

    def test_failing_destination_does_not_stop_others(self):
        """ A destination with a live run is skipped, and one failing still lets the rest run """
        synced = []

        def fake_run_sync(sync_destination):
            synced.append(sync_destination)
            if sync_destination == self.failing:
                raise ConnectionError("remote down")

        with patch("sync.tasks.sync_tasks.run_sync", side_effect=fake_run_sync), \
                patch("sync.tasks.sync_tasks.report_exc_info"):
            with self.assertRaises(RuntimeError):
                sync_all()
        self.assertCountEqual(synced, [self.failing, self.working])
