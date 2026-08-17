"""
#1701: the recovery command must find exactly the runs whose VEP output went missing, and leave
everything else finished.
"""
from io import StringIO

from django.core.management import call_command
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationRangeLock, AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

VEP_WARNINGS_ONE_SKIP = ("WARNING: line 3 skipped (1 200 . GG A . . variant_id=1...): "
                         "Length of reference allele (GG length 2) does not match coordinates 200-200")


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class FixTruncatedAnnotationRunsTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)

    def _finished_run(self, **kwargs) -> AnnotationRun:
        """ A run that reached FINISHED - upload_end with no error_exception is what get_status() needs """
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)
        kwargs.setdefault("pipeline_type", VariantAnnotationPipelineType.STANDARD)
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=lock,
                                                      upload_end=timezone.now(), **kwargs)
        self.assertEqual(annotation_run.status, AnnotationStatus.FINISHED)
        return annotation_run

    def _run_command(self, *args) -> str:
        out = StringIO()
        call_command("fix_truncated_annotation_runs", *args, stdout=out)
        return out.getvalue()

    def test_reports_only_affected_runs(self):
        truncated = self._finished_run(dump_count=12501, annotated_count=0)
        complete = self._finished_run(dump_count=500, annotated_count=500)
        skips_accounted_for = self._finished_run(dump_count=500, annotated_count=499,
                                                 vep_warnings=VEP_WARNINGS_ONE_SKIP)
        not_counted = self._finished_run()  # No dump_count/annotated_count - nothing to check

        output = self._run_command()
        reported = {int(line.split("\t")[0]) for line in output.splitlines()
                    if line.split("\t")[0].isdigit()}
        self.assertEqual(reported, {truncated.pk})
        for annotation_run in (complete, skips_accounted_for, not_counted):
            self.assertNotIn(annotation_run.pk, reported)
        self.assertIn("1 affected runs, 12501 variants missing annotation", output)

    def test_min_shortfall_floor(self):
        self._finished_run(dump_count=500, annotated_count=499)  # Could just be a CSQ-less record
        big = self._finished_run(dump_count=12501, annotated_count=0)

        output = self._run_command("--min-shortfall", "10")
        reported = {int(line.split("\t")[0]) for line in output.splitlines()
                    if line.split("\t")[0].isdigit()}
        self.assertEqual(reported, {big.pk})

    def test_mark_error_moves_only_affected_runs(self):
        truncated = self._finished_run(dump_count=12501, annotated_count=0)
        complete = self._finished_run(dump_count=500, annotated_count=500)

        self._run_command("--mark-error")

        truncated.refresh_from_db()
        complete.refresh_from_db()
        self.assertEqual(truncated.status, AnnotationStatus.ERROR)
        self.assertIn("#1701", truncated.error_exception)
        self.assertIn("shortfall 12501", truncated.error_exception)
        self.assertEqual(complete.status, AnnotationStatus.FINISHED)

    def test_report_does_not_change_status(self):
        truncated = self._finished_run(dump_count=12501, annotated_count=0)
        self._run_command()
        truncated.refresh_from_db()
        self.assertEqual(truncated.status, AnnotationStatus.FINISHED)
