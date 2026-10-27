"""download_annotation_data: an optional index (.csi) is fetched only when variantgrid.com hosts it (#2499), a
plugin we ship ("refresh") is re-fetched when the installed copy differs from the hosted one, a file that ships with
the repo is never downloaded, and --bash-script covers the same files as a download unless --download-all."""
from io import StringIO
from unittest import mock

from django.core.management import call_command
from django.test import SimpleTestCase

_COMMAND = "variantopedia.management.commands.download_annotation_data"
_PLUGIN_URL = "https://raw.githubusercontent.com/SACGF/variantgrid/refs/heads/master/x/UTRAnnotator.pm"
_CHECKS = {
    "GRCh38/gnomad4": {"valid": True, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz"},
    "GRCh38/gnomad4_tbi": {"valid": False, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz.tbi"},
    "GRCh38/gnomad4_csi": {"valid": False, "fix": "cd /data/x;wget https://variantgrid.com/download/x/g.vcf.gz.csi",
                           "optional": True},
    "GRCh38/plugin_UTRAnnotator": {"valid": True, "fix": f"cd /data/plugins;wget {_PLUGIN_URL}", "refresh": True},
    "GRCh38/cytoband": {"valid": False, "fix": "git checkout -- variantgrid/data/reference/hg38/cytoband.hg38.txt.gz"},
}


class DownloadAnnotationDataOptionalIndexTest(SimpleTestCase):
    def _dry_run(self, hosted: bool, plugin_differs=False) -> str:
        out = StringIO()
        with mock.patch(f"{_COMMAND}.annotation_data_exists", return_value=_CHECKS), \
                mock.patch(f"{_COMMAND}._is_hosted", return_value=hosted) as is_hosted, \
                mock.patch(f"{_COMMAND}._differs_from_hosted", return_value=plugin_differs) as differs:
            call_command("download_annotation_data", "--dry-run", stdout=out)
        # only the optional file is probed; a required file is always listed as missing
        is_hosted.assert_called_once_with("https://variantgrid.com/download/x/g.vcf.gz.csi")
        # only the plugin we ship is compared with its hosted copy
        differs.assert_called_once_with(_PLUGIN_URL, "/data/plugins/UTRAnnotator.pm")
        return out.getvalue()

    def test_shipped_plugin_is_refetched_only_when_it_differs(self):
        self.assertNotIn("plugin_UTRAnnotator", self._dry_run(hosted=True))
        output = self._dry_run(hosted=True, plugin_differs=True)
        self.assertIn("Found 3 missing or outdated", output)
        self.assertIn("GRCh38/plugin_UTRAnnotator", output)

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

    def test_repo_shipped_file_is_never_downloaded(self):
        output = self._dry_run(hosted=True)
        self.assertNotIn("cytoband", output)
        self.assertNotIn("Could not parse", output)

    def _bash_script(self, *args, plugin_differs=False) -> list[str]:
        out = StringIO()
        with mock.patch(f"{_COMMAND}.annotation_data_exists", return_value=_CHECKS), \
                mock.patch(f"{_COMMAND}._differs_from_hosted", return_value=plugin_differs):
            call_command("download_annotation_data", "--bash-script", *args, stdout=out)
        script = out.getvalue()
        self.assertNotIn("cytoband", script)
        return [line for line in script.splitlines() if line.startswith("wget")]

    def test_bash_script_downloads_only_missing_and_outdated_files(self):
        lines = self._bash_script()
        self.assertEqual(len(lines), 2)
        self.assertNotIn("||", lines[0])
        self.assertTrue(lines[1].endswith("|| echo 'optional, not hosted: GRCh38/gnomad4_csi'"))
        # a shipped plugin overwrites the installed copy
        lines = self._bash_script(plugin_differs=True)
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[2], f"wget -O /data/plugins/UTRAnnotator.pm {_PLUGIN_URL}")

    def test_bash_script_download_all_covers_every_hosted_file(self):
        lines = self._bash_script("--download-all")
        self.assertEqual(len(lines), 4)
        self.assertTrue(lines[0].endswith("https://variantgrid.com/download/x/g.vcf.gz"))
        self.assertEqual(lines[3], f"wget -O /data/plugins/UTRAnnotator.pm {_PLUGIN_URL}")
