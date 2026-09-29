"""
Tests for the get_partition_archive_path helper.

@see claude/issue_1537_archive_plan.md §1
"""

import os
import tempfile

from django.conf import settings
from django.test import TestCase, override_settings

from library.django_utils.partition_archive_paths import get_partition_archive_path


class PartitionArchivePathTests(TestCase):

    def test_path_under_db_subdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            with override_settings(PARTITION_ARCHIVE_DIR=tmp):
                path = get_partition_archive_path("foo.dump")

            db_name = settings.DATABASES["default"]["NAME"]
            self.assertEqual(path, os.path.join(tmp, db_name, "foo.dump"))
            self.assertTrue(os.path.isdir(os.path.join(tmp, db_name)))
