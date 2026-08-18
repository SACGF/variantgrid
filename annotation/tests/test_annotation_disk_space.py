"""
Tests for annotation disk-space handling (issue #1670).

Two halves:
 - cleanup: every route into annotation.signals.annotation_run_cleanup, which is the only place a run's
   on-disk output is removed. A failed run keeps everything so it can be investigated.
 - gating: the dispatcher declines to start new VEP runs when the dump filesystem is nearly full, but
   keeps importing runs already past VEP - which is what frees the space again.
"""
import os
import tempfile
from datetime import timedelta
from unittest import mock

from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.signals.annotation_run_cleanup import remove_annotation_run_output
from annotation.signals.manual_signals import annotation_run_complete_signal
from annotation.tasks import annotation_scheduler_task
from annotation.annotation_run_files import get_annotsv_dir
from annotation.pipelines import get_runner
from annotation.pipelines.vep import VEPRunner
from annotation.tasks.annotate_variants import (
    _cleanup_reclaimed_run_files,
    annotate_variants,
    annotation_run_retry,
    import_annotation_run,
)
from annotation.tasks.annotation_scheduler_task import (
    count_annotation_run,
    dispatch_annotation_runs,
    has_free_disk_for_annotation,
)
from genes.models_enums import AnnotationConsortium
from library.utils.file_utils import DiskUsage
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

STANDARD = VariantAnnotationPipelineType.STANDARD


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationRunCleanupTestCase(TestCase):
    """ Each way into the cleanup module, and the one case that must not reach it. """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def _run_with_output(self, tmp_dir, pipeline_type=STANDARD, **run_kwargs):
        """ A run whose full set of pipeline artifacts exists on disk, as its own runner names them.
            Returns (run, paths, annotsv_dir) - annotsv_dir is None for pipelines that don't use one. """
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)
        run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type)
        run.vcf_dump_filename = os.path.join(tmp_dir, f"dump_{run.pk}__task.vcf")
        runner = get_runner(pipeline_type)

        annotsv_dir = None
        if pipeline_type == VariantAnnotationPipelineType.ANNOTSV:
            annotsv_dir = get_annotsv_dir(run)
            os.makedirs(annotsv_dir, exist_ok=True)
            # AnnotSV writes more than just the TSV into its directory
            open(os.path.join(annotsv_dir, "AnnotSV.log"), "w").close()

        paths = runner.get_output_paths(run, run.vcf_dump_filename)
        for path in paths:
            open(path, "w").close()
        runner.record_resume_state(run, run.vcf_dump_filename)
        for k, v in run_kwargs.items():
            setattr(run, k, v)
        run.save()
        return run, paths, annotsv_dir

    def _assert_all_removed(self, paths, annotsv_dir):
        for path in paths:
            self.assertFalse(os.path.exists(path), f"left {os.path.basename(path)} behind")
        if annotsv_dir:
            self.assertFalse(os.path.exists(annotsv_dir))

    def test_success_signal_reclaims_output(self):
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir,
                                  ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS=True):
            run, paths, annotsv_dir = self._run_with_output(
                tmp_dir, pipeline_type=VariantAnnotationPipelineType.ANNOTSV)

            annotation_run_complete_signal.send(sender="test", annotation_run=run,
                                                variant_annotation_version=self.vav,
                                                pipeline_type=run.pipeline_type)

            self._assert_all_removed(paths, annotsv_dir)

    def test_success_cleanup_can_be_turned_off(self):
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir,
                                  ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS=False):
            run, paths, annotsv_dir = self._run_with_output(
                tmp_dir, pipeline_type=VariantAnnotationPipelineType.ANNOTSV)

            annotation_run_complete_signal.send(sender="test", annotation_run=run,
                                                variant_annotation_version=self.vav,
                                                pipeline_type=run.pipeline_type)

            for path in paths:
                self.assertTrue(os.path.exists(path))
            self.assertTrue(os.path.isdir(annotsv_dir))

    def test_success_signal_without_a_run_is_harmless(self):
        # The signal predates the cleanup - upload's receiver sends/consumes it without a run, so the
        # handler's None guard is load-bearing. No assertions: this passes if nothing raises.
        annotation_run_complete_signal.send(sender="test", variant_annotation_version=self.vav,
                                            pipeline_type=STANDARD)

    def test_failed_import_keeps_everything(self):
        # The point of #1670: cleanup hangs off the success signal, so a failure has no route into the
        # cleanup module at all and every file survives for investigation.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir,
                                  ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS=True):
            run, paths, _ = self._run_with_output(tmp_dir)

            with mock.patch.object(VEPRunner, "import_results",
                            side_effect=RuntimeError("import blew up")), \
                    mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
                with self.assertRaises(RuntimeError):
                    import_annotation_run.apply((run.pk,)).get()

            for path in paths:
                self.assertTrue(os.path.exists(path), f"failure removed {os.path.basename(path)}")

    def test_discarded_attempt_keeps_shared_annotsv_dir(self):
        # #1658: the AnnotSV dir is keyed on the run and shared between attempts, so an attempt that lost
        # the run must not take the winner's output with it.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run, paths, annotsv_dir = self._run_with_output(
                tmp_dir, pipeline_type=VariantAnnotationPipelineType.ANNOTSV)

            _cleanup_reclaimed_run_files(run)

            for path in paths:
                self.assertFalse(os.path.exists(path))
            self.assertTrue(os.path.isdir(annotsv_dir))
            self.assertTrue(os.path.exists(os.path.join(annotsv_dir, "AnnotSV.log")))

    def test_queryset_delete_removes_output_files(self):
        # post_delete rather than a delete() override, so a bulk delete (eg
        # subdivide_annotation_range_lock) cannot silently orphan the files.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run, paths, annotsv_dir = self._run_with_output(
                tmp_dir, pipeline_type=VariantAnnotationPipelineType.ANNOTSV)

            with mock.patch.object(AnnotationRun, "delete_related_objects"):
                AnnotationRun.objects.filter(pk=run.pk).delete()

            self._assert_all_removed(paths, annotsv_dir)

    def test_retry_upload_only_rejects_run_whose_vcf_was_cleaned_up(self):
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run, _, _ = self._run_with_output(tmp_dir)
            remove_annotation_run_output(run)

            with self.assertRaises(ValueError) as cm:
                annotation_run_retry(run, upload_only=True)
            self.assertIn("no longer exists", str(cm.exception))


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationDiskGateTestCase(TestCase):
    """ The dispatcher's low-disk brake. """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 200000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(4)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def _make_run(self, lo_idx, hi_idx, status=AnnotationStatus.CREATED):
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[lo_idx],
                                                  max_variant=self.variants[hi_idx], count=100)
        run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STANDARD)
        if status == AnnotationStatus.ANNOTATION_COMPLETED:
            now = timezone.now()  # past VEP, waiting on the import lane
            run.dump_start = now - timedelta(minutes=5)
            run.dump_end = now - timedelta(minutes=4)
            run.dump_count = 100
            run.annotation_start = now - timedelta(minutes=4)
            run.annotation_end = now
            run.vcf_annotated_filename = "/does/not/need/to/exist.vcf.gz"
            run.save()
            self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        return run

    @staticmethod
    def _disk_usage(free_gigs):
        return DiskUsage(mount_point="/data", available_kb=int(free_gigs * 1000000),
                              available_nice=f"{free_gigs}G", percent_nice="95%")

    def _dispatch(self, free_gigs):
        with mock.patch.object(annotation_scheduler_task, "annotation_worker_slots", return_value=4), \
                mock.patch.object(annotation_scheduler_task, "_annotation_dump_disk_usage",
                                  return_value=self._disk_usage(free_gigs)), \
                mock.patch.object(annotate_variants, "apply_async") as vep_launch, \
                mock.patch.object(import_annotation_run, "apply_async") as import_launch, \
                mock.patch.object(count_annotation_run, "apply_async"):
            dispatch_annotation_runs(self.vav.pk)
        return vep_launch, import_launch

    @override_settings(ANNOTATION_MIN_FREE_DISK_GIGS=10)
    def test_low_disk_blocks_vep_but_still_imports(self):
        self._make_run(0, 1)  # CREATED - needs a dump + VEP, so needs disk
        self._make_run(2, 3, status=AnnotationStatus.ANNOTATION_COMPLETED)  # past VEP - frees disk

        vep_launch, import_launch = self._dispatch(free_gigs=2.0)

        vep_launch.assert_not_called()
        self.assertEqual(import_launch.call_count, 1)

    @override_settings(ANNOTATION_MIN_FREE_DISK_GIGS=10)
    def test_ample_disk_dispatches_vep(self):
        self._make_run(0, 1)
        vep_launch, _ = self._dispatch(free_gigs=500.0)
        self.assertEqual(vep_launch.call_count, 1)

    @override_settings(ANNOTATION_MIN_FREE_DISK_GIGS=None)
    def test_gate_disabled_dispatches_regardless(self):
        self._make_run(0, 1)
        vep_launch, _ = self._dispatch(free_gigs=0.1)
        self.assertEqual(vep_launch.call_count, 1)

    @override_settings(ANNOTATION_MIN_FREE_DISK_GIGS=10)
    def test_unmeasurable_disk_does_not_block(self):
        # Not being able to measure is a reason to keep annotating and let warn_low_disk_space fire,
        # not to silently halt the pipeline.
        with mock.patch.object(annotation_scheduler_task, "_annotation_dump_disk_usage", return_value=None):
            self.assertTrue(has_free_disk_for_annotation())
