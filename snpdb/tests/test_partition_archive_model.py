"""
Tests for PartitionArchive model + archive_partitioned_model helper.

@see claude/issue_1537_archive_plan.md §7
"""

import tempfile
from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotationVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.models.models_partition_archive import PartitionArchive
from snpdb.partition_archive import (
    PartitionArchivePreconditionError,
    archive_partitioned_model,
)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class PartitionArchiveHelperTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.user = User.objects.create_user(username="archiver")

    def _new_vav(self):
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        return VariantAnnotationVersion.objects.create(
            **kwargs, status=VariantAnnotationVersion.Status.HISTORICAL
        )

    def _archive_dir_override(self):
        return override_settings(PARTITION_ARCHIVE_DIR=self._tmpdir.name)

    def setUp(self):
        super().setUp()
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)

    def test_resolve_source_returns_live_row_or_none(self):
        from annotation.models import AnnotationVersion
        vav = self._new_vav()
        with self._archive_dir_override():
            with mock.patch("snpdb.partition_archive._schedule_archive_task"):
                archive = archive_partitioned_model(vav, self.user)
        self.assertIsNotNone(archive.resolve_source())
        vav_pk = vav.pk
        AnnotationVersion.objects.filter(variant_annotation_version_id=vav_pk).delete()
        VariantAnnotationVersion.objects.filter(pk=vav_pk).delete()
        self.assertIsNone(archive.resolve_source())

    def test_rejects_model_without_data_archive_mixin(self):
        with self._archive_dir_override():
            with self.assertRaises(PartitionArchivePreconditionError):
                archive_partitioned_model(self.user, self.user)

    def test_rejects_already_archived_source(self):
        vav = self._new_vav()
        vav.data_archived_date = timezone.now()
        vav.save()
        with self._archive_dir_override():
            with self.assertRaises(PartitionArchivePreconditionError):
                archive_partitioned_model(vav, self.user)

    def test_rejects_concurrent_in_progress_archive(self):
        vav = self._new_vav()
        with self._archive_dir_override():
            with mock.patch("snpdb.partition_archive._schedule_archive_task"):
                archive_partitioned_model(vav, self.user)
            with self.assertRaises(PartitionArchivePreconditionError):
                archive_partitioned_model(vav, self.user)

    def test_successful_schedule_creates_pending_row(self):
        vav = self._new_vav()
        with self._archive_dir_override():
            with mock.patch("snpdb.partition_archive._schedule_archive_task") as sched:
                with self.captureOnCommitCallbacks(execute=True):
                    archive = archive_partitioned_model(vav, self.user, reason="why")
            sched.assert_called_once_with(archive.pk, "why")

        self.assertEqual(archive.status, PartitionArchive.Status.PENDING)
        self.assertEqual(archive.source_app_label, "annotation")
        self.assertEqual(archive.source_model, "VariantAnnotationVersion")
        self.assertEqual(archive.source_pk, str(vav.pk))
        self.assertIn("annotation_variantannotation",
                      archive.source_table_names[0])
        self.assertEqual(archive.dumped_by, self.user)
