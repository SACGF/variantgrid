import csv
import os
import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from analysis.enums import TagLocation
from analysis.models import VariantTag
from annotation.fake_data import create_fake_variant_annotation, get_fake_annotation_version
from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from library.django_utils.django_partition import temporary_db_table
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.enums import ImportSource
from snpdb.models import GenomeBuild, Tag, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_test_variant,
    slowly_create_test_variant_from_coordinate,
)
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
        self._import_csv([TAGS_GRID_CSV_HEADER,
                          ["3:128198980 G>A", "GRCh37", "GATA2", tag_id, "", "", "", "",
                           self.user.username, "2026-10-29 02:46:19.957469+00:00", self.variant.pk, 1]])

    def _import_csv(self, rows: list[list]):
        tmp_dir = tempfile.mkdtemp()
        csv_filename = os.path.join(tmp_dir, "tags.csv")
        with open(csv_filename, "w", newline="") as f:
            csv.writer(f).writerows(rows)

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

    def test_tags_grid_download_round_trips_symbolic_variants(self, _mock_liftover):
        """ An annotated variant re-imports from its g.HGVS. Without annotation it falls back to the Variant
            string, which writes a symbolic variant with its end, the form VariantCoordinate.from_string reads """
        unannotated_deletion = slowly_create_test_variant_from_coordinate(
            VariantCoordinate(chrom="21", position=47532744, ref="T", alt="<DEL>", svlen=-1224), self.grch37)
        annotated_deletion = slowly_create_test_variant_from_coordinate(
            VariantCoordinate(chrom="3", position=128200125, ref="G", alt="<DEL>", svlen=-1201), self.grch37)
        vav = get_fake_annotation_version(self.grch37).variant_annotation_version
        partition_table = vav.get_partition_table(base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition_table):
            create_fake_variant_annotation(annotated_deletion, vav)

        exported_tag = Tag.objects.create(pk="exported")
        for variant in (self.variant, unannotated_deletion, annotated_deletion):
            variant_tag = VariantTag.objects.create(variant=variant, tag=exported_tag, genome_build=self.grch37,
                                                    user=self.user)
            assign_permission_to_user_and_groups(self.user, variant_tag)

        self.client.force_login(self.user)
        url = reverse("variant_tags_datatable", kwargs={"genome_build_name": self.grch37.name})
        response = self.client.get(url, {"tag": exported_tag.pk, "dataTableCsv": 1})
        header, *rows = csv.reader(b"".join(response.streaming_content).decode().splitlines())
        variant_col = header.index("Variant")
        hgvs_g_col = header.index("g.HGVS")
        rows_by_hgvs_g = {row[hgvs_g_col]: row for row in rows}
        self.assertIn("21:47532744-47533968 <DEL>", [row[variant_col] for row in rows])
        # Unreadable, so only the g.HGVS can bring this one back
        rows_by_hgvs_g["NC_000003.11:g.128200130_128201330del"][variant_col] = "unreadable"

        self._import_csv([header, *rows])
        imported = VariantTag.objects.filter(tag=exported_tag, location=TagLocation.EXTERNAL)
        self.assertEqual(set(imported.values_list("variant", flat=True)),
                         {self.variant.pk, unannotated_deletion.pk, annotated_deletion.pk})
