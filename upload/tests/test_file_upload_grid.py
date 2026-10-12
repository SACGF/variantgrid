from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase

from snpdb.models.models_enums import ImportSource
from snpdb.views.datatable_view import datatable_response
from upload.grids import FileUploadColumns
from upload.models import FileUpload, UploadedFileTypes


class FileUploadColumnsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="file_upload_grid_user")
        for name, file_type in [("mine.vcf", UploadedFileTypes.VCF),
                                ("liftover.vcf", UploadedFileTypes.LIFTOVER),
                                ("genes.txt", UploadedFileTypes.GENE_LIST)]:
            FileUpload.objects.create(user=cls.user, name=name, file_type=file_type,
                                      import_source=ImportSource.WEB_UPLOAD)

    def _names(self, **params) -> set[str]:
        request = RequestFactory().get("/", params)
        request.user = self.user
        return {row["name"]["text"] for row in datatable_response(FileUploadColumns(request))["data"]}

    def test_internal_types_hidden_by_default(self):
        self.assertEqual(self._names(), {"mine.vcf", "genes.txt"})

    def test_file_type_select(self):
        self.assertEqual(self._names(file_type=UploadedFileTypes.LIFTOVER), {"liftover.vcf"})
        self.assertEqual(self._names(file_type=FileUploadColumns.ALL_FILE_TYPES),
                         {"mine.vcf", "liftover.vcf", "genes.txt"})
