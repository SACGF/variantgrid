"""annotation_files_check: a file under ANNOTATION_BASE_DIR is fetched from variantgrid.com; one that ships with the
repo (the cytobands) is restored with git, never wget."""
from django.test import SimpleTestCase, override_settings

from variantgrid.deployment_validation.annotation_files_check import _get_fix_instructions


@override_settings(ANNOTATION_BASE_DIR="/data/annotation", BASE_DIR="/opt/variantgrid")
class FixInstructionsTest(SimpleTestCase):
    def test_hosted_file_is_a_wget(self):
        self.assertEqual(_get_fix_instructions("/data/annotation/fasta/GRCh38.fna.gz"),
                         "cd /data/annotation/fasta;wget https://variantgrid.com/download/annotation/fasta/GRCh38.fna.gz")

    def test_repo_file_is_a_checkout(self):
        fix = _get_fix_instructions("/opt/variantgrid/variantgrid/data/reference/hg38/cytoband.hg38.txt.gz")
        self.assertEqual(fix, "git checkout -- variantgrid/data/reference/hg38/cytoband.hg38.txt.gz")
