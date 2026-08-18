from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models.models_enums import ImportSource
from upload.models import FileUpload, UploadedFileTypes


class TestUploadedFilePermissions(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="file_owner", password="x")
        cls.other_user = User.objects.create_user(username="file_other", password="x")
        cls.superuser = User.objects.create_superuser(username="file_super", password="x")

        cls.file_upload = FileUpload.objects.create(
            user=cls.owner,
            name="test_file.vcf",
            path="/tmp/test_file.vcf",
            file_type=UploadedFileTypes.VCF,
            import_source=ImportSource.COMMAND_LINE,
        )

    def test_can_view_is_owner_or_superuser(self):
        """ FileUpload doesn't use Guardian - it's owner-or-superuser only """
        self.assertTrue(self.file_upload.can_view(self.owner))
        self.assertTrue(self.file_upload.can_view(self.superuser))
        self.assertFalse(self.file_upload.can_view(self.other_user))

