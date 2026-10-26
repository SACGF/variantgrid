import os
import tempfile

from django.test import SimpleTestCase, override_settings

from variantgrid.deployment_validation.unused_annotation_files import (
    UnusedAnnotationData,
    unused_annotation_data,
)


class UnusedAnnotationDataTest(SimpleTestCase):

    def setUp(self):
        tmp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(tmp_dir.cleanup)
        self.vep_dir = tmp_dir.name
        self.data_dir = os.path.join(self.vep_dir, "annotation_data")

    def _file(self, rel_path: str, size: int) -> str:
        path = os.path.join(self.data_dir, rel_path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as f:
            f.write(b"x" * size)
        return path

    def test_reports_what_the_enabled_builds_do_not_name(self):
        self._file("GRCh37/dbNSFP5.gz", 10)
        self._file("GRCh37/dbNSFP5.gz.tbi", 1)
        old_dbnsfp = self._file("GRCh37/dbNSFP4.gz", 20)
        self._file("GRCh37/dbNSFP4.gz.tbi", 2)
        self._file("GRCh37/af_greater_than_5.vcf.gz", 3)
        old_dir = os.path.dirname(self._file("GRCh37/old/gnomad.vcf.gz", 5))
        self._file("GRCh37/old/gnomad.vcf.gz.tbi", 1)
        self._file("all_builds/maxentscan/splice5sequences", 4)
        self._file("GRCh38/dbNSFP4.gz", 100)  # GRCh38 isn't enabled, so its directory isn't looked at

        annotation = {
            "GRCh37": {"enabled": True, "columns_version": 5, "vep_config": {
                "sift": True,
                "dbnsfp": "annotation_data/GRCh37/dbNSFP5.gz",
                "maxentscan": "annotation_data/all_builds/maxentscan",
                "cosmic": None,
            }},
            "GRCh38": {"enabled": False, "vep_config": {"dbnsfp": "annotation_data/GRCh38/dbNSFP5.gz"}},
        }
        common_filters = {"GRCh37": {"gnomad_af_filename": "annotation_data/GRCh37/af_greater_than_5.vcf.gz"}}
        with override_settings(ANNOTATION_VEP_BASE_DIR=self.vep_dir, ANNOTATION=annotation,
                               VCF_IMPORT_COMMON_FILTERS=common_filters):
            unused = unused_annotation_data()

        # The index counts with its data file, and a directory holding nothing used is one entry
        self.assertEqual({os.path.join(self.data_dir, "GRCh37"): [UnusedAnnotationData(old_dbnsfp, 22),
                                                                   UnusedAnnotationData(old_dir, 6)]}, unused)
