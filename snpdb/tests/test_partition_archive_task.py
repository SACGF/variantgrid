"""
Tests for the perform_partition_archive Celery task and admin lifecycle helpers.

@see claude/issue_1537_archive_plan.md §7
"""

import os
import subprocess
import tempfile
from unittest import mock

from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase, override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotationVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.models.models_partition_archive import PartitionArchive
from snpdb.partition_archive import (
    archive_partitioned_model,
    clear_partition_archive_dump,
    mark_partition_archive_restored,
)
from snpdb.tasks.partition_archive_tasks import perform_partition_archive


def _fake_subprocess_run_factory(dump_path: str, table_listing: list[str], dump_should_fail: bool = False,
                                 listing_missing: bool = False):
    """ Return a side_effect for subprocess.run that fakes pg_dump and pg_restore --list. """

    def _run(cmd, *args, **kwargs):
        binary = cmd[0]
        if binary == "pg_dump":
            if dump_should_fail:
                raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="boom")
            with open(dump_path, "wb") as f:
                f.write(b"PGDMP-fake")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
        if binary == "pg_restore":
            present = [] if listing_missing else table_listing
            return subprocess.CompletedProcess(
                cmd, 0,
                stdout="\n".join(f"; TABLE - public {t}" for t in present),
                stderr="",
            )
        raise AssertionError(f"unexpected subprocess invocation: {cmd}")

    return _run


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class PerformPartitionArchiveTaskTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.user = User.objects.create_user(username="archiver")

    def setUp(self):
        super().setUp()
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.archive_dir_override = override_settings(PARTITION_ARCHIVE_DIR=self._tmpdir.name)
        self.archive_dir_override.enable()
        self.addCleanup(self.archive_dir_override.disable)

        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        self.vav = VariantAnnotationVersion.objects.create(
            **kwargs, status=VariantAnnotationVersion.Status.HISTORICAL
        )
        self.child_tables = [
            self.vav.get_partition_table(base_table_name=t)
            for t in self.vav.RECORDS_BASE_TABLE_NAMES
        ]
        with mock.patch("snpdb.partition_archive._schedule_archive_task"):
            self.archive = archive_partitioned_model(self.vav, self.user, reason="r")

    def _patch_subprocess(self, **kwargs):
        return mock.patch(
            "snpdb.tasks.partition_archive_tasks.subprocess.run",
            side_effect=_fake_subprocess_run_factory(self.archive.dump_path,
                                                    self.child_tables, **kwargs),
        )

    def test_happy_path(self):
        with self._patch_subprocess():
            perform_partition_archive(self.archive.pk, "r")

        self.archive.refresh_from_db()
        self.assertEqual(self.archive.status, PartitionArchive.Status.COMPLETE)
        self.assertIsNotNone(self.archive.sha256)
        self.assertIsNotNone(self.archive.completed_at)

        existing_tables = set(connection.introspection.table_names())
        for child in self.child_tables:
            self.assertNotIn(child, existing_tables)

        self.vav.refresh_from_db()
        self.assertIsNotNone(self.vav.data_archived_date)
        self.assertEqual(self.vav.data_archived_by_id, self.user.pk)
        self.assertEqual(self.vav.data_restorable_from, self.archive.dump_path)

        from eventlog.models import Event
        self.assertTrue(Event.objects.filter(name="partition_archive_complete").exists())

    def test_pg_dump_failure(self):
        with self._patch_subprocess(dump_should_fail=True):
            perform_partition_archive(self.archive.pk, "r")

        self.archive.refresh_from_db()
        self.assertEqual(self.archive.status, PartitionArchive.Status.FAILED)
        self.assertTrue(self.archive.error_message)

        existing_tables = set(connection.introspection.table_names())
        for child in self.child_tables:
            self.assertIn(child, existing_tables)

        self.vav.refresh_from_db()
        self.assertIsNone(self.vav.data_archived_date)

    def test_pg_restore_listing_missing_table(self):
        with self._patch_subprocess(listing_missing=True):
            perform_partition_archive(self.archive.pk, "r")

        self.archive.refresh_from_db()
        self.assertEqual(self.archive.status, PartitionArchive.Status.FAILED)

        existing_tables = set(connection.introspection.table_names())
        for child in self.child_tables:
            self.assertIn(child, existing_tables)

    def test_skips_when_already_in_progress(self):
        self.archive.status = PartitionArchive.Status.IN_PROGRESS
        self.archive.save()
        with mock.patch("snpdb.tasks.partition_archive_tasks.subprocess.run") as run:
            perform_partition_archive(self.archive.pk, "r")
        run.assert_not_called()


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class PartitionArchiveAdminHelperTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="archiver")

    def _make_complete_archive(self, vav=None) -> PartitionArchive:
        if vav is None:
            grch37 = GenomeBuild.get_name_or_alias("GRCh37")
            kwargs = get_fake_vep_version(grch37, AnnotationConsortium.ENSEMBL, 2)
            vav = VariantAnnotationVersion.objects.create(
                **kwargs, status=VariantAnnotationVersion.Status.HISTORICAL,
            )
            vav.data_archived_date = timezone.now()
            vav.data_archived_by = self.user
            vav.data_archive_reason = "t"
            vav.data_restorable_from = "/tmp/x.dump"
            vav.save()
        archive = PartitionArchive.objects.create(
            archive_name=f"annotation_variantannotationversion_{vav.pk}_x",
            dump_path="/tmp/x.dump",
            source_app_label="annotation",
            source_model="VariantAnnotationVersion",
            source_pk=str(vav.pk),
            source_table_names=["annotation_variantannotation"],
            status=PartitionArchive.Status.COMPLETE,
            dumped_by=self.user,
        )
        return archive

    def test_mark_restored_rejects_non_complete(self):
        archive = self._make_complete_archive()
        archive.status = PartitionArchive.Status.PENDING
        archive.save()
        with self.assertRaises(ValueError):
            mark_partition_archive_restored(archive, self.user)

    def test_mark_restored_clears_data_archive_fields(self):
        archive = self._make_complete_archive()
        mark_partition_archive_restored(archive, self.user)
        archive.refresh_from_db()
        self.assertIsNotNone(archive.restored_date)
        vav = VariantAnnotationVersion.objects.get(pk=int(archive.source_pk))
        self.assertIsNone(vav.data_archived_date)
        self.assertIsNone(vav.data_restorable_from)

    def test_clear_dump_deletes_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(b"x")
            dump_path = f.name
        archive = self._make_complete_archive()
        archive.dump_path = dump_path
        archive.save()

        clear_partition_archive_dump(archive, self.user)
        archive.refresh_from_db()
        self.assertIsNotNone(archive.cleared_date)
        self.assertFalse(os.path.exists(dump_path))
