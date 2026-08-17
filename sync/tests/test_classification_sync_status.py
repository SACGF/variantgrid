from django.contrib.auth.models import User
from django.test import Client, TestCase, override_settings
from django.urls import reverse

from classification.enums import ShareLevel, SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.models import Lab
from sync.classification_sync_status import ClassificationSyncState, classification_sync_status
from sync.models.enums import SyncStatus
from sync.models.models import SyncDestination, SyncRun
from sync.models.models_classification_sync import ClassificationModificationSyncRecord
from sync.shariant.variant_grid_upload import VariantGridUploadSyncer

REMOTE_HOST = "https://shariant.org.au"
SYNC_DETAILS = {"test_shariant": {"host": REMOTE_HOST}}


@override_settings(CLASSIFICATION_MATCH_VARIANTS=False,
                   SYNC_DETAILS=SYNC_DETAILS,
                   STATICFILES_STORAGE='django.contrib.staticfiles.storage.StaticFilesStorage')
class ClassificationSyncStatusTestCase(TestCase):

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
                "filters": {"allele_origin": ["germline"]}
            }
        )

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def _create_classification(
            self,
            share_level: ShareLevel = ShareLevel.ALL_USERS,
            allele_origin: str = "germline",
            lab: Lab = None,
            user: User = None) -> Classification:
        vc = Classification.create(
            user=user or self.user,
            lab=lab or self.lab,
            lab_record_id=None,
            data={
                "allele_origin": allele_origin,
                "clinical_significance": "VUS"
            },
            save=True,
            source=SubmissionSource.API
        )
        vc.publish_latest(user or self.user, share_level=share_level)
        return vc

    def _create_sync_record(self, vc: Classification, remote_pk: int = 999):
        sync_run = SyncRun.objects.create(destination=self.sync_destination, status=SyncStatus.SUCCESS)
        return ClassificationModificationSyncRecord.objects.create(
            run=sync_run,
            classification_modification=vc.last_published_version,
            meta={"meta": {"id": remote_pk}}
        )

    def _only_status(self, vc: Classification):
        statuses = classification_sync_status(vc)
        self.assertEqual(len(statuses), 1, "Expected one status for the one upload destination")
        return statuses[0]

    def test_synced(self):
        vc = self._create_classification()
        self._create_sync_record(vc, remote_pk=1234)

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.SYNCED)
        self.assertEqual(status.remote_url, f"{REMOTE_HOST}/classification/classification/1234")
        self.assertIsNotNone(status.last_synced)
        self.assertFalse(status.reasons)

    def test_changes_pending(self):
        vc = self._create_classification()
        self._create_sync_record(vc, remote_pk=1234)

        vc.patch_value({"clinical_significance": "LP"}, user=self.user, save=True, source=SubmissionSource.API)
        vc.publish_latest(self.user, share_level=ShareLevel.ALL_USERS)

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.CHANGES_PENDING)
        self.assertEqual(status.remote_url, f"{REMOTE_HOST}/classification/classification/1234")

    def test_pending_falls_back_to_remote_search(self):
        """ Remotes yet to be upgraded have no view_classification_lab_record to link to """
        vc = self._create_classification()

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.PENDING)
        self.assertEqual(status.remote_url, f"{REMOTE_HOST}/variantopedia/search?search={vc.lab_record_id}")
        self.assertIsNone(status.last_synced)

    def test_pending_links_to_lab_record_once_remote_is_upgraded(self):
        self.sync_destination.config["remote_lab_record_url"] = True
        self.sync_destination.save()
        vc = self._create_classification()

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.PENDING)
        self.assertEqual(
            status.remote_url,
            f"{REMOTE_HOST}/classification/classification/lab_record/instx/labby/{vc.lab_record_id}")

    def test_lab_share_level_has_no_statuses(self):
        """ The share level message already says the record stays with the lab """
        vc = self._create_classification(share_level=ShareLevel.LAB)

        self.assertEqual(classification_sync_status(vc), [])

    def test_share_level_is_still_an_exclusion_reason(self):
        vc = self._create_classification(share_level=ShareLevel.LAB)

        runner = VariantGridUploadSyncer()
        runner.configure(self.sync_destination)
        reasons = runner.exclusion_reasons(vc.last_published_version)
        self.assertTrue(any("uploads once shared more widely" in reason for reason in reasons),
                        f"Expected a share level reason in {reasons}")

    def test_excluded_by_lab_mapping(self):
        self.sync_destination.config["mapping"]["labs"] = {"instx/some_other_lab": True}
        self.sync_destination.save()
        vc = self._create_classification()

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertTrue(any("is not configured to share with" in reason for reason in status.reasons),
                        f"Expected a lab mapping reason in {status.reasons}")
        self.assertIsNone(status.remote_url)

    def test_excluded_by_filter_clause(self):
        vc = self._create_classification(allele_origin="somatic")

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertIn("Allele origin is Somatic, shares when set to Germline", status.reasons)

    def test_excluded_by_somatic_filter_clause(self):
        self.sync_destination.config["filters"] = {"somatic": False}
        self.sync_destination.save()
        vc = self._create_classification(allele_origin="somatic")

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertIn("Allele origin is Somatic, shares unless set to Somatic", status.reasons)

    def test_blank_value_reported_for_filter_clause(self):
        vc = self._create_classification(allele_origin=None)

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertIn("Allele origin is blank, shares when set to Germline", status.reasons)

    def test_excluded_after_an_earlier_upload_notes_the_remote_copy_is_stale(self):
        vc = self._create_classification()
        self._create_sync_record(vc, remote_pk=1234)

        vc.patch_value({"allele_origin": "somatic"}, user=self.user, save=True, source=SubmissionSource.API)
        vc.publish_latest(self.user, share_level=ShareLevel.ALL_USERS)

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertEqual(status.remote_url, f"{REMOTE_HOST}/classification/classification/1234")
        self.assertIn("would not be sent", status.note)

    def test_no_stale_remote_note_when_nothing_was_uploaded(self):
        vc = self._create_classification(allele_origin="somatic")

        status = self._only_status(vc)
        self.assertEqual(status.state, ClassificationSyncState.EXCLUDED)
        self.assertIsNone(status.note)

    def test_excluded_by_filter_clause_with_label_override(self):
        self.sync_destination.config["filter_labels"] = {"allele_origin": "Germline records upload to Shariant"}
        self.sync_destination.save()
        vc = self._create_classification(allele_origin="somatic")

        status = self._only_status(vc)
        self.assertEqual(status.reasons, ["Germline records upload to Shariant"])

    def test_external_lab_has_no_statuses(self):
        ext_lab = Lab.objects.get(group_name='instx/ext')
        ext_user = User.objects.get(username='joejoe2')
        vc = self._create_classification(lab=ext_lab, user=ext_user)

        self.assertEqual(classification_sync_status(vc), [])

    def test_unpublished_record_has_no_statuses(self):
        vc = Classification.create(
            user=self.user, lab=self.lab, lab_record_id=None,
            data={"allele_origin": "germline"}, save=True, source=SubmissionSource.API)

        self.assertEqual(classification_sync_status(vc), [])

    def test_disabled_destination_has_no_statuses(self):
        self.sync_destination.enabled = False
        self.sync_destination.save()
        vc = self._create_classification()

        self.assertEqual(classification_sync_status(vc), [])

    def test_classification_page_renders_the_sharing_card(self):
        vc = self._create_classification(allele_origin="somatic")

        client = Client()
        client.force_login(self.user)
        response = client.get(reverse('view_classification', kwargs={'classification_id': vc.pk}))

        self.assertEqual(response.status_code, 200)
        content = response.content.decode()
        self.assertIn("Test Shariant", content)
        self.assertIn("Allele origin is Somatic, shares when set to Germline", content)

    def test_reporting_leaves_filter_config_alone(self):
        """ convert_to_q_w_key strips nulls out of list values, which would quietly change what syncs """
        self.sync_destination.config["filters"] = {"allele_origin": ["germline", None]}
        self.sync_destination.save()
        vc = self._create_classification(allele_origin="somatic")

        runner = VariantGridUploadSyncer()
        runner.configure(self.sync_destination)
        runner.exclusion_reasons(vc.last_published_version)
        self.assertEqual(runner.filters["allele_origin"], ["germline", None])
