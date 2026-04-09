from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models.models_enums import ImportSource
from upload.models import UploadedFile, UploadedFileTypes


class TestUploadedFilePermissions(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="file_owner", password="x")
        cls.other_user = User.objects.create_user(username="file_other", password="x")
        cls.superuser = User.objects.create_superuser(username="file_super", password="x")

        cls.uploaded_file = UploadedFile.objects.create(
            user=cls.owner,
            name="test_file.vcf",
            path="/tmp/test_file.vcf",
            file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.COMMAND_LINE,
        )

    def test_owner_can_view(self):
        self.assertTrue(self.uploaded_file.can_view(self.owner))

    def test_non_owner_cannot_view(self):
        self.assertFalse(self.uploaded_file.can_view(self.other_user))

    def test_superuser_can_view_any_file(self):
        self.assertTrue(self.uploaded_file.can_view(self.superuser))

