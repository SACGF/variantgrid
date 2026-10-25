from unittest.mock import MagicMock, patch

from django.test import TestCase

from classification.tests.models.test_utils import ClassificationTestUtils
from library.oauth import ServerAuth
from snpdb.models.models import Lab
from sync.shariant.variant_grid_download import shariant_download_to_upload


@patch("sync.shariant.variant_grid_download.report_message")
class ShariantDownloadToUploadTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.server_auth = ServerAuth(username="user", password="password", host="https://shariant.example.com")
        self.sync_destination = MagicMock()
        self.sync_destination.name = "Test Shariant"

    def _map(self, record_json: dict):
        return shariant_download_to_upload(self.server_auth, self.sync_destination, MagicMock(), record_json)

    def test_own_lab_skipped_by_id_not_meta(self, _report_message):
        """ meta.lab_id names an external lab, but the id (which the inserter files under) is our own lab's """
        record = {"id": "instx/labby/vc1", "meta": {"id": 1, "lab_id": "instx/ext"}, "data": {}}
        self.assertIsNone(self._map(record))

    def test_own_lab_delete_skipped(self, _report_message):
        self.assertIsNone(self._map({"id": "instx/labby/vc1", "delete": True}))

    def test_external_lab_delete_passed_on(self, _report_message):
        upload = self._map({"id": "instx/ext/vc1", "delete": True})
        self.assertEqual(upload["id"], "instx/ext/vc1")
        self.assertTrue(upload["delete"])

    def test_unknown_lab(self, _report_message):
        self.assertIsNone(self._map({"id": "neworg/newlab/vc1", "delete": True}))
        self.assertFalse(Lab.objects.filter(group_name="neworg/newlab").exists())

        upload = self._map({"id": "neworg/newlab/vc1", "meta": {"id": 1, "lab_name": "New Lab"}, "data": {}})
        self.assertEqual(upload["id"], "neworg/newlab/vc1")
        self.assertTrue(Lab.objects.get(group_name="neworg/newlab").external)

    def test_id_without_lab_skipped(self, report_message):
        for record_id in ["1234", "1234.5678", None, 1234, "bad org/lab/vc1"]:
            self.assertIsNone(self._map({"id": record_id, "delete": True}), msg=repr(record_id))
        self.assertEqual(report_message.call_count, 5)
