import csv
import os
import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from analysis.enums import TagLocation
from analysis.models import VariantTag
from snpdb.enums import ImportSource
from snpdb.models import GenomeBuild, Tag
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from upload.enums import UploadedFileTypes
from upload.models import FileUpload, UploadPipeline, UploadStep
from upload.tasks.import_variant_tags_task import VariantTagsCreateVCFTask, VariantTagsInsertTask

# Header written by the Tags grid CSV download (variantopedia.grids.VariantTagsColumns)
TAGS_GRID_CSV_HEADER = ["Variant", "Genome Build", "Gene", "Tag", "Analysis ID", "Analysis", "Sample", "Patient",
                        "Username", "Created", "Variant ID", "Tag ID"]


@patch("upload.tasks.import_variant_tags_task.create_liftover_pipelines")
class ImportVariantTagsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="import_variant_tags_user", password="x")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variant = slowly_create_test_variant("3", 128198980, "G", "A", cls.grch37)
        # An existing tag imports whatever its name - only new tag ids must be alphanumeric
        Tag.objects.create(pk="existing tag")

    def _import(self, tag_id: str):
        tmp_dir = tempfile.mkdtemp()
        csv_filename = os.path.join(tmp_dir, "tags.csv")
        with open(csv_filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(TAGS_GRID_CSV_HEADER)
            writer.writerow(["3:128198980 G>A", "GRCh37", "GATA2", tag_id, "", "", "", "",
                             self.user.username, "2026-10-29 02:46:19.957469+00:00", self.variant.pk, 1])

        file_upload = FileUpload.objects.create(user=self.user, name="tags.csv", path=csv_filename,
                                                file_type=UploadedFileTypes.VARIANT_TAGS,
                                                import_source=ImportSource.COMMAND_LINE)
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        create_step = UploadStep.objects.create(upload_pipeline=upload_pipeline, name="create", sort_order=0,
                                                input_filename=csv_filename,
                                                output_filename=os.path.join(tmp_dir, "tags.vcf"))
        VariantTagsCreateVCFTask.process_items(create_step)
        insert_step = UploadStep.objects.create(upload_pipeline=upload_pipeline, name="insert", sort_order=1)
        VariantTagsInsertTask.process_items(insert_step)

    def test_tags_grid_download_imports(self, _mock_liftover):
        self._import("existing tag")
        variant_tag = VariantTag.objects.get(tag_id="existing tag")
        self.assertEqual(variant_tag.variant, self.variant)
        self.assertEqual(variant_tag.location, TagLocation.EXTERNAL)
        self.assertEqual(variant_tag.user, self.user)

    def test_new_tag_must_be_alphanumeric(self, _mock_liftover):
        with self.assertRaisesRegex(ValueError, "must be alphanumeric"):
            self._import("new tag")
        self.assertFalse(Tag.objects.filter(pk="new tag").exists())
