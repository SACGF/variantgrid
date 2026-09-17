"""
The one-off sweep of directories that leaked before #928 - manage.py import_processing_cleanup.

The rule it has to get right is "only remove what its owner is finished with", so each case here is a
directory whose owner is present-and-done, present-and-still-working, or gone entirely.
"""
import os
import tempfile
import time
from io import StringIO

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase, override_settings

from annotation.models.models import ManualVariantEntryCollection
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from upload.models import (
    FileUpload,
    ProcessingStatus,
    UploadedFileTypes,
    UploadedManualVariantEntryCollection,
    UploadPipeline,
)

OLD_MTIME = time.time() - 30 * 24 * 60 * 60


class ImportProcessingCleanupTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="import_processing_cleanup_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _make_pipeline(self, status: str) -> UploadPipeline:
        mvec = ManualVariantEntryCollection.objects.create(user=self.user, genome_build=self.grch37)
        file_upload = FileUpload.objects.create(user=self.user, name="manual_variant_entry",
                                                path="/nowhere/manual_variant_entry.vcf",
                                                file_type=UploadedFileTypes.MANUAL_VARIANT_ENTRY,
                                                import_source=ImportSource.WEB)
        UploadedManualVariantEntryCollection.objects.create(file_upload=file_upload, collection=mvec)
        return UploadPipeline.objects.create(status=status, file_upload=file_upload)

    @staticmethod
    def _make_dir(import_processing_dir: str, name: str, empty=False) -> str:
        path = os.path.join(import_processing_dir, name)
        os.makedirs(path)
        if not empty:
            with open(os.path.join(path, "scratch.txt"), "w") as f:
                f.write("x")
        os.utime(path, (OLD_MTIME, OLD_MTIME))
        return path

    def test_removes_only_what_owners_are_done_with(self):
        successful = self._make_pipeline(ProcessingStatus.SUCCESS)
        errored = self._make_pipeline(ProcessingStatus.ERROR)
        with tempfile.TemporaryDirectory() as import_processing_dir:
            paths = {
                "successful": self._make_dir(import_processing_dir, f"pipeline_{successful.pk}"),
                "errored": self._make_dir(import_processing_dir, f"pipeline_{errored.pk}"),
                "no_row": self._make_dir(import_processing_dir, "pipeline_99999999"),
                "clingen_empty": self._make_dir(import_processing_dir, "clingen_allele_registry_abc", empty=True),
                "clingen_dump": self._make_dir(import_processing_dir, "clingen_allele_registry_failures"),
                "unit_test": self._make_dir(import_processing_dir, "test"),
                "unknown": self._make_dir(import_processing_dir, "something_we_dont_know"),
            }
            recent = os.path.join(import_processing_dir, "pipeline_99999998")
            os.makedirs(recent)

            with override_settings(IMPORT_PROCESSING_DIR=import_processing_dir):
                call_command("import_processing_cleanup", dry_run=True, stdout=StringIO())
                self.assertTrue(all(os.path.exists(p) for p in paths.values()), "--dry-run removed something")

                call_command("import_processing_cleanup", stdout=StringIO())

            for key in ("successful", "no_row", "clingen_empty", "unit_test"):
                self.assertFalse(os.path.exists(paths[key]), f"{key} should have been removed")
            for key in ("errored", "clingen_dump", "unknown"):
                self.assertTrue(os.path.exists(paths[key]), f"{key} should have been kept")
            self.assertTrue(os.path.exists(recent), "an entry modified today should have been kept")
