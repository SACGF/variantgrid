from django.test import TestCase

from library.utils.file_utils import get_extension_without_gzip


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
