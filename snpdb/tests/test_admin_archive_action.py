"""
Tests for the admin "Archive partition data" action.

@see claude/issue_1537_archive_plan.md §7
"""

import tempfile
from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotationVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.models.models_partition_archive import PartitionArchive


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AdminArchiveActionTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.admin = User.objects.create_superuser(
            username="admin", email="a@example.com", password="x",
        )

    def _new_vav(self):
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, 2)
        return VariantAnnotationVersion.objects.create(
            **kwargs, status=VariantAnnotationVersion.Status.HISTORICAL
        )

    def test_admin_action_creates_pending_archive(self):
        vav = self._new_vav()
        self.client.force_login(self.admin)
        url = reverse("admin:annotation_variantannotationversion_changelist")

        with tempfile.TemporaryDirectory() as tmp:
            with override_settings(PARTITION_ARCHIVE_DIR=tmp):
                with mock.patch("snpdb.partition_archive._schedule_archive_task") as sched:
                    with self.captureOnCommitCallbacks(execute=True):
                        response = self.client.post(url, {
                            "action": "archive_partition_data",
                            "_selected_action": [vav.pk],
                        }, follow=True)

        self.assertEqual(response.status_code, 200)
        archive = PartitionArchive.objects.get(source_pk=str(vav.pk))
        self.assertEqual(archive.status, PartitionArchive.Status.PENDING)
        sched.assert_called_once_with(archive.pk, "Admin action")
