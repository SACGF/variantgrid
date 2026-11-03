from unittest.mock import patch

from django.test import TestCase
from django.test.utils import override_settings

from annotation.enums import VariantAnnotationPipelineType
from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.vep_annotation import get_vep_command
from annotation.vep_command_check import check_vep_command, normalise_vep_command, vep_command_diffs
from snpdb.models import GenomeBuild

STANDARD = VariantAnnotationPipelineType.STANDARD


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class VEPCommandCheckTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.vav = get_fake_annotation_version(cls.grch37).variant_annotation_version

    def _cmd(self, vcf="in.vcf", output="out.vcf"):
        return get_vep_command(vcf, output, self.grch37, self.vav.annotation_consortium, STANDARD,
                               variant_annotation_version=self.vav)

    def test_normalised_command_ignores_run_and_box(self):
        command = normalise_vep_command(self._cmd())
        with override_settings(ANNOTATION_VEP_MEMORY_LIMIT_GB=4, ANNOTATION_VEP_FORK=8,
                               ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT="/path/perlbrew_runner.sh"):
            other_run = normalise_vep_command(self._cmd("/dump/run_2.vcf", "/dump/run_2.vep.vcf.gz"))
        self.assertEqual(command, other_run)
        self.assertTrue(command.startswith("vep "))

    @patch("annotation.vep_command_check.report_message")
    def test_first_run_recorded_later_change_reported_not_raised(self, mock_report_message):
        check_vep_command(self.vav, STANDARD, self._cmd())
        self.vav.refresh_from_db()
        recorded = self.vav.vep_commands[STANDARD]
        mock_report_message.assert_not_called()

        with override_settings(ANNOTATION_VEP_PICK_ORDER="mane_select,canonical"):
            check_vep_command(self.vav, STANDARD, self._cmd())
            self.assertEqual(1, len(vep_command_diffs(self.vav)))
        mock_report_message.assert_called_once()
        self.vav.refresh_from_db()
        self.assertEqual(recorded, self.vav.vep_commands[STANDARD], "Keeps what the version first ran")
