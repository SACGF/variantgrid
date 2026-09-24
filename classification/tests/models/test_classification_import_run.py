from unittest.mock import patch

from django.test import TestCase
from django.urls import reverse

from classification.enums import SubmissionSource
from classification.models.classification import Classification
from classification.models.classification_import_run import ClassificationImportRun
from classification.models.classification_inserter import BulkClassificationInserter
from classification.models.uploaded_classifications_unmapped import (
    UploadedClassificationsUnmapped,
    UploadedClassificationsUnmappedStatus,
)
from classification.tests.models.test_utils import ClassificationTestUtils


class ClassificationImportRunAlreadyWithdrawnTest(TestCase):

    def setUp(self):
        ClassificationTestUtils.setUp()

    @patch("classification.models.classification_import_run.MAX_LAB_RECORD_IDS_ALREADY_WITHDRAWN", 1)
    def test_lab_record_ids_already_withdrawn(self):
        lab, user = ClassificationTestUtils.lab_and_user()
        for lab_record_id in ["withdrawn_1", "withdrawn_2"]:
            classification = Classification.create(user=user, lab=lab, lab_record_id=lab_record_id,
                                                   data={"condition": "x"}, save=True,
                                                   source=SubmissionSource.API)
            classification.publish_latest(user=user)
            classification.set_withdrawn(user=user, withdraw=True)

        upload = UploadedClassificationsUnmapped.objects.create(
            url="file:///test.json", filename="test.json", lab=lab, user=user,
            status=UploadedClassificationsUnmappedStatus.Processed, validation_summary={"message_counts": {}})
        import_run = ClassificationImportRun.objects.create(identifier="test_already_withdrawn", from_file=upload)
        inserter = BulkClassificationInserter(user=user)
        for lab_record_id in ["withdrawn_1", "withdrawn_2"]:
            response = inserter.insert({"id": f"{lab.group_name}/{lab_record_id}", "upsert": {"condition": "y"}},
                                       submission_source=SubmissionSource.API, import_run=import_run)
            import_run.increment_status(response)
        inserter.finish()
        import_run.save()

        self.assertEqual(import_run.row_count_already_withdrawn, 2)
        self.assertEqual(import_run.lab_record_ids_already_withdrawn, ["withdrawn_1"], "Capped at max")
        self.assertTrue(import_run.lab_record_ids_already_withdrawn_truncated)

        self.client.force_login(user)
        response = self.client.get(reverse("classification_upload_unmapped_status_detail",
                                           kwargs={"uploaded_classification_unmapped_id": upload.pk}))
        withdrawn_1 = Classification.objects.get(lab=lab, lab_record_id="withdrawn_1")
        self.assertContains(response, reverse("view_classification", kwargs={"classification_id": withdrawn_1.pk}))
        self.assertContains(response, "Only the first 1 are listed")
