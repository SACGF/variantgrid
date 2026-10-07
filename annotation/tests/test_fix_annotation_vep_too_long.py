""" One-off backfill of TOO_LONG rows (#2104) - delete along with fix_annotation_vep_too_long """
from annotation.management.commands.fix_annotation_vep_too_long import fix_annotation_vep_too_long
from annotation.models.models_enums import VEPSkippedReason
from annotation.tests.test_vep_too_long import VEPTooLongTestBase


class FixAnnotationVEPTooLongTestCase(VEPTooLongTestBase):
    def test_backfill_finished_run(self):
        sv_run = self._sv_run(self.long_sv)
        sv_run.dump_count = 0
        sv_run.save()

        fix_annotation_vep_too_long()
        self.assertEqual(self._skipped_reasons(sv_run), {self.long_sv.pk: VEPSkippedReason.TOO_LONG})
