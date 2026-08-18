"""
#720: AnnotSV as its own annotation pipeline type.

Covers what the split buys: a supplementary run waits for the VEP run whose rows it updates, does not
hold up a VCF import, sees SVs too large for VEP, and appears for every existing range lock the moment
the pipeline is enabled - which is the backfill the issue asked for.
"""
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_vep_version
from annotation.models import AnnotationVersion, VariantAnnotationVersion
from annotation.models.models import AnnotationRangeLock, AnnotationRun
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.pipelines import blocking_pipeline_types, enabled_pipeline_types, get_runner
from annotation.tasks.annotation_scheduler_task import (
    _dispatchable_runs_qs,
    _handle_variant_annotation_version,
    count_annotation_run,
)
from genes.models_enums import AnnotationConsortium
from library.utils import execute_cmd
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

ANNOTSV = VariantAnnotationPipelineType.ANNOTSV
STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


class AnnotSVPipelineTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def _make_lock(self):
        return AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)

    def _make_runs(self, lock):
        sv_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STRUCTURAL)
        annotsv_run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=ANNOTSV)
        return sv_run, annotsv_run

    def test_annotsv_run_waits_for_sv_vep_run(self):
        lock = self._make_lock()
        sv_run, annotsv_run = self._make_runs(lock)
        dispatchable = _dispatchable_runs_qs(self.vav, timezone.now())
        self.assertIn(sv_run.pk, set(dispatchable.values_list("pk", flat=True)))
        self.assertNotIn(annotsv_run.pk, set(dispatchable.values_list("pk", flat=True)))

        sv_run.dump_count = 0
        sv_run.save()
        self.assertEqual(sv_run.status, AnnotationStatus.FINISHED)
        dispatchable = _dispatchable_runs_qs(self.vav, timezone.now())
        self.assertIn(annotsv_run.pk, set(dispatchable.values_list("pk", flat=True)))

    def test_sv_free_range_finishes_annotsv_run_without_dispatching(self):
        """ What makes enabling AnnotSV cheap: most range locks hold no SVs at all, so the count lane
            finishes their AnnotSV runs on db_workers without one ever reaching a worker. Counting is
            deliberately not gated on the dependency - it asks how many SVs are in range, which does not
            depend on VEP having run. """
        lock = self._make_lock()  # only the two SNVs from setUpTestData
        _sv_run, annotsv_run = self._make_runs(lock)

        count_annotation_run(annotsv_run.pk)

        annotsv_run.refresh_from_db()
        self.assertEqual(annotsv_run.count, 0)
        self.assertEqual(annotsv_run.status, AnnotationStatus.FINISHED)

    def test_annotsv_does_not_block_vcf_import(self):
        self.assertNotIn(ANNOTSV, blocking_pipeline_types())
        self.assertIn(STANDARD, blocking_pipeline_types())
        self.assertIn(STRUCTURAL, blocking_pipeline_types())

    @override_settings(ANNOTATION_ANNOTSV_ENABLED=True)
    def test_enabling_annotsv_creates_runs_for_existing_locks(self):
        """ #720's backfill: enabling the pipeline makes the scheduler's orphaned-lock repair create an
            ANNOTSV run for every lock that already finished VEP. """
        lock = self._make_lock()
        for pipeline_type in (STANDARD, STRUCTURAL, VariantAnnotationPipelineType.GENE_LEVEL):
            run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type)
            run.dump_count = 0
            run.save()
        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).exists())

        _handle_variant_annotation_version(self.vav)

        self.assertTrue(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                     pipeline_type=ANNOTSV).exists())

    def test_annotsv_runs_not_created_while_disabled(self):
        self.assertNotIn(ANNOTSV, enabled_pipeline_types())
        lock = self._make_lock()
        _handle_variant_annotation_version(self.vav)
        self.assertFalse(AnnotationRun.objects.filter(annotation_range_lock=lock,
                                                      pipeline_type=ANNOTSV).exists())

    @override_settings(ANNOTATION_VEP_SV_MAX_SIZE=1000)
    def test_annotsv_does_not_inherit_veps_sv_size_cap(self):
        """ AnnotSV was fed the VEP dump, so it never saw SVs above ANNOTATION_VEP_SV_MAX_SIZE - the ones
            its ACMG ranking is for. Its own dump applies no such cap. """
        lock = self._make_lock()
        sv_run, annotsv_run = self._make_runs(lock)
        # VEPRunner adds the abs_svlen annotation solely to apply the cap
        self.assertIn("abs_svlen", get_runner(STRUCTURAL).get_variants_qs(sv_run).query.annotations)
        self.assertNotIn("abs_svlen", get_runner(ANNOTSV).get_variants_qs(annotsv_run).query.annotations)


class ExecuteCmdTimeoutTest(TestCase):
    """ #720: AnnotSV needs a timeout on execute_cmd before it can move off subprocess.run and onto the
        lease heartbeat's process_callback. """

    def test_timeout_kills_and_returns_non_zero(self):
        return_code, _std_out, _std_err = execute_cmd(["sleep", "30"], timeout=0.2)
        self.assertNotEqual(return_code, 0)

    def test_no_timeout_returns_normally(self):
        return_code, std_out, _std_err = execute_cmd(["echo", "hello"], timeout=30)
        self.assertEqual(return_code, 0)
        self.assertEqual((std_out or "").strip(), "hello")
