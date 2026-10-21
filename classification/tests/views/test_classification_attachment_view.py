from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from classification.enums import SubmissionSource
from classification.models.classification import Classification
from classification.models.variant_models import ClassificationAttachment
from classification.tests.models.test_utils import ClassificationTestUtils
from library.enums.file_attachments import AttachmentFileType


class ClassificationAttachmentViewTestCase(TestCase):
    """ An attachment is only served to users who can view its classification """

    def setUp(self):
        ClassificationTestUtils.setUp()
        lab, self.lab_user = ClassificationTestUtils.lab_and_user()
        self.outsider = User.objects.get_or_create(username="attachment_outsider")[0]
        classification = Classification.create(user=self.lab_user, lab=lab, lab_record_id="attachment_test",
                                               data={}, save=True, source=SubmissionSource.API)
        attachment = ClassificationAttachment.objects.create(
            classification=classification, file_type=AttachmentFileType.TEXT,
            file=SimpleUploadedFile("attachment_test.txt", b"lab only"))
        self.addCleanup(attachment.file.delete, save=False)
        self.url = reverse("view_classification_file_attachment", kwargs={"pk": attachment.pk})
        self.thumbnail_url = reverse("view_classification_file_attachment_thumbnail", kwargs={"pk": attachment.pk})

    def tearDown(self):
        ClassificationTestUtils.tearDown()

    def test_outsider_is_refused_file_and_thumbnail(self):
        self.client.force_login(self.outsider)
        self.assertEqual(self.client.get(self.url).status_code, 403)
        self.assertEqual(self.client.get(self.thumbnail_url).status_code, 403)

    def test_lab_user_gets_file(self):
        self.client.force_login(self.lab_user)
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content, b"lab only")
