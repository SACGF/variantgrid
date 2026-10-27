"""download_annotation_data: an optional index (.csi) is fetched only when variantgrid.com hosts it (#2499)."""
from io import StringIO
from unittest import mock

from django.core.management import call_command
from django.test import SimpleTestCase

_CHECKS = {
    "GRCh38/gnomad4": {"valid": True, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz"},
    "GRCh38/gnomad4_tbi": {"valid": False, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz.tbi"},
    "GRCh38/gnomad4_csi": {"valid": False, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz.csi",
                           "optional": True},
}


class DownloadAnnotationDataOptionalIndexTest(SimpleTestCase):
    def _dry_run(self, hosted: bool) -> str:
        out = StringIO()
        with mock.patch("variantopedia.management.commands.download_annotation_data.annotation_data_exists",
                        return_value=_CHECKS), \
                mock.patch("variantopedia.management.commands.download_annotation_data._is_hosted",
                           return_value=hosted) as is_hosted:
            call_command("download_annotation_data", "--dry-run", stdout=out)
        # only the optional file is probed; a required file is always listed as missing
        is_hosted.assert_called_once_with("https://variantgrid.com/download/x/g.vcf.gz.csi")
        return out.getvalue()

    def test_hosted_csi_is_listed_with_the_missing_files(self):
        output = self._dry_run(hosted=True)
        self.assertIn("Found 2 missing", output)
        self.assertIn("GRCh38/gnomad4_csi", output)
        self.assertNotIn("not hosted", output)

    def test_unhosted_csi_is_skipped_not_failed(self):
        output = self._dry_run(hosted=False)
        self.assertIn("Found 1 missing", output)
        self.assertIn("GRCh38/gnomad4_tbi", output)
        self.assertIn("not hosted on variantgrid.com, skipped", output)
        self.assertIn("build_tabix_csi.sh", output)

    def test_bash_script_tolerates_a_missing_optional_file(self):
        out = StringIO()
        with mock.patch("variantopedia.management.commands.download_annotation_data.annotation_data_exists",
                        return_value=_CHECKS):
            call_command("download_annotation_data", "--bash-script", stdout=out)
        lines = [line for line in out.getvalue().splitlines() if line.startswith("wget")]
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[2].endswith("|| echo 'optional, not hosted: GRCh38/gnomad4_csi'"))
        self.assertNotIn("||", lines[1])
