from django.test import SimpleTestCase, TestCase

from library.utils.file_utils import get_extension_without_gzip, get_mount_point_for_directory


class TestGetExtensionWithoutGzip(TestCase):
    """Tests for get_extension_without_gzip — drives all file-type routing in the upload pipeline."""

    def test_vcf_gz(self):
        self.assertEqual(get_extension_without_gzip("sample.vcf.gz"), "vcf")

    def test_vcf_bgz(self):
        """bgz is a separate extension used by htslib/tabix-indexed files."""
        self.assertEqual(get_extension_without_gzip("sample.vcf.bgz"), "vcf")

    def test_plain_vcf(self):
        self.assertEqual(get_extension_without_gzip("sample.vcf"), "vcf")

    def test_no_extension_returns_empty_string(self):
        self.assertEqual(get_extension_without_gzip("myfile"), "")


class MountPointMatchingTestCase(SimpleTestCase):
    """ Which filesystem a directory is attributed to. """

    MOUNTS = ["/", "/data", "/database", "/boot/efi"]

    def test_longest_match_wins(self):
        # Every absolute path is under '/', so a plain prefix test would attribute this to whichever df
        # listed first.
        self.assertEqual(get_mount_point_for_directory("/data/annotation_scratch", self.MOUNTS), "/data")

    def test_separator_boundary_respected(self):
        # '/database'.startswith('/data') is True, but they are different filesystems
        self.assertEqual(get_mount_point_for_directory("/database/dumps", self.MOUNTS), "/database")

    def test_mount_point_itself(self):
        self.assertEqual(get_mount_point_for_directory("/data", self.MOUNTS), "/data")

    def test_falls_back_to_root(self):
        self.assertEqual(get_mount_point_for_directory("/home/user/variantgrid", self.MOUNTS), "/")

    def test_no_match(self):
        self.assertIsNone(get_mount_point_for_directory("/data/x", ["/mnt"]))
