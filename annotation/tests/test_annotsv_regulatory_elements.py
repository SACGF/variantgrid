import os
import shutil
import subprocess
import tempfile

from django.conf import settings
from django.test import SimpleTestCase

from annotation.annotsv_regulatory_elements import is_bed_sorted, sort_command


class AnnotSVRegulatoryElementsTest(SimpleTestCase):
    """ #2749: the fixture is broken the way the Annotations_Human 3.5 ABC file is - starts ascend within
        every block, but chromosome 2 comes back after its alt contig, which bedtools -sorted rejects. """
    FIXTURE = os.path.join(settings.BASE_DIR, "annotation/tests/test_data/annotsv",
                           "regulatory_elements_chrom_reappears.sorted.bed")

    def test_reappearing_chromosome_is_unsorted(self):
        self.assertFalse(is_bed_sorted(self.FIXTURE))

    def test_printed_command_sorts_the_file_in_place(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            # A directory an operator's shell has to quote
            filename = os.path.join(tmp_dir, "regulatory elements.sorted.bed")
            shutil.copy(self.FIXTURE, filename)

            subprocess.run(sort_command(filename), shell=True, check=True)

            self.assertTrue(is_bed_sorted(filename))
            with open(self.FIXTURE, encoding="utf-8") as before, open(filename, encoding="utf-8") as after:
                self.assertEqual(sorted(before), sorted(after))
            self.assertEqual(os.listdir(tmp_dir), ["regulatory elements.sorted.bed"])
