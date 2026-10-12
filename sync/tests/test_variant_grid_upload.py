from unittest.mock import MagicMock, patch

from django.test import TestCase, override_settings

from classification.enums import ShareLevel, SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from sync.models.models import SyncDestination
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.variant_grid_upload import VariantGridUploadSyncer
from sync.sync_runner import SyncRunInstance

SYNC_DETAILS = {"test_shariant": {"host": "https://shariant.org.au"}}


@override_settings(CLASSIFICATION_MATCH_VARIANTS=False, SYNC_DETAILS=SYNC_DETAILS)
class VariantGridUploadSyncerTestCase(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.sync_destination = SyncDestination.objects.create(
            name="Test Shariant",
            config={
                "type": "shariant",
                "direction": "upload",
                "sync_details": "test_shariant",
                "mapping": {"labs": {"instx/labby": True}},
            }
        )
        self.classification = Classification.create(
            user=self.user,
            lab=self.lab,
            lab_record_id=None,
            data={"allele_origin": "germline", "clinical_significance": "VUS"},
            save=True,
            source=SubmissionSource.API
        )
        self.classification.publish_latest(self.user, share_level=ShareLevel.ALL_USERS)

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def _delta_records(self) -> list:
        syncer = VariantGridUploadSyncer()
        syncer.configure(self.sync_destination)
        return list(syncer.records_to_sync())

    def _sync(self, remote_result: dict) -> list[dict]:
        """ Runs a delta upload with the remote answering remote_result for every record, returns the records sent """
        response = MagicMock()
        response.json.side_effect = lambda: {"results": [remote_result] * len(sent_records)}
        server_auth = MagicMock()
        sent_records = []

        def post(**kwargs):
            sent_records.extend(kwargs["json"]["records"])
            return response

        server_auth.post.side_effect = post
        with patch.object(SyncRunInstance, "server_auth", return_value=server_auth):
            VariantGridUploadSyncer().sync(SyncRunInstance(self.sync_destination))
        return sent_records

    def test_internal_error_is_not_synced(self):
        self._sync({"internal_error": "boom"})

        sync_record = ClassificationModificationSyncRecord.objects.get()
        self.assertFalse(sync_record.success)
        self.assertEqual(self._delta_records(), [self.classification.last_published_version])

    def test_withdrawal_sent_on_delta(self):
        self._sync({"meta": {"id": 1234}})
        self.assertEqual(self._delta_records(), [])

        self.classification.set_withdrawn(self.user, withdraw=True)
        sent_records = self._sync({"meta": {"id": 1234}, "withdrawn": True})
        self.assertEqual(len(sent_records), 1)
        self.assertTrue(sent_records[0].get("delete"))

        # the remote has answered for the withdrawal, so it isn't sent again
        self.assertEqual(self._delta_records(), [])
