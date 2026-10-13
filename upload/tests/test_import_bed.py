"""A BED file whose genome build can't be worked out waits for the user to set it (#2088)"""
import os
import tempfile
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase, override_settings

from snpdb.forms import GenomicIntervalsCollectionForm
from snpdb.models.models_enums import ImportSource, ImportStatus, ProcessingStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.views.datatable_view import datatable_response
from upload.grids import FileUploadColumns
from upload.models import FileUpload, UploadedBed, UploadedFileTypes, UploadPipeline
from upload.tasks.import_bedfile_task import ImportBedFileTask
from upload.upload_processing import process_uploaded_file

NO_HEADER_BED = os.path.join(settings.BASE_DIR, "upload", "test_data", "bed", "test_no_header.bed")


@patch.object(ImportBedFileTask, "_get_genome_build", return_value=None)
class TestImportBedRequiresGenomeBuild(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="bed_import_user", password="x")

    def setUp(self):
        super().setUp()
        processed_bed_dir = tempfile.TemporaryDirectory()
        self.addCleanup(processed_bed_dir.cleanup)
        settings_override = override_settings(PROCESSED_BED_FILES_DIR=processed_bed_dir.name)
        settings_override.enable()
        self.addCleanup(settings_override.disable)

    def _import_bed(self) -> UploadPipeline:
        file_upload = FileUpload.objects.create(user=self.user, name="test_no_header.bed", path=NO_HEADER_BED,
                                                file_type=UploadedFileTypes.BED,
                                                import_source=ImportSource.COMMAND_LINE)
        upload_pipeline, _ = process_uploaded_file(file_upload, run_async=False)
        return UploadPipeline.objects.get(pk=upload_pipeline.pk)

    def test_pipeline_waits_for_user_input(self, _get_genome_build):
        upload_pipeline = self._import_bed()
        self.assertEqual(ProcessingStatus.TERMINATED_EARLY, upload_pipeline.status, upload_pipeline.progress_status)

        gic = UploadedBed.objects.get(file_upload=upload_pipeline.file_upload).genomic_intervals_collection
        self.assertEqual(ImportStatus.REQUIRES_USER_INPUT, gic.import_status)

        request = RequestFactory().get("/")
        request.user = self.user
        [row] = datatable_response(FileUploadColumns(request))["data"]
        self.assertEqual(gic.get_absolute_url(), row["status"]["url"])

    @patch("upload.models.models_uploaded_files.process_bed_file", return_value=2)  # bedtools isn't on CI
    def test_setting_genome_build_finishes_pipeline(self, _process_bed_file, _get_genome_build):
        upload_pipeline = self._import_bed()
        gic = UploadedBed.objects.get(file_upload=upload_pipeline.file_upload).genomic_intervals_collection

        data = {"name": gic.name, "user": self.user.pk, "genome_build": GenomeBuild.grch37().pk}
        form = GenomicIntervalsCollectionForm(data, instance=gic)
        self.assertTrue(form.is_valid(), form.errors)
        form.save()

        gic.refresh_from_db()
        self.assertEqual(ImportStatus.SUCCESS, gic.import_status)
        upload_pipeline.refresh_from_db()
        self.assertEqual(ProcessingStatus.SUCCESS, upload_pipeline.status)
        self.assertEqual("Success", upload_pipeline.progress_status)
        self.assertEqual(2, upload_pipeline.items_processed)
