"""
Tests for the capacity-limited, merging annotation dispatcher (issue #2667).

Mirrors analysis/tests/test_scheduler.py: the scheduler now only creates pending AnnotationRuns and
dispatch_annotation_runs is the single authority that leases + launches them, merging the backlog
into bigger batches when (and only when) workers are saturated.

annotate_variants.apply_async and annotation_worker_slots are mocked so we assert dispatch decisions
without running VEP or inspecting a live celery pool.
"""
from contextlib import ExitStack
from datetime import timedelta
from unittest import mock

from celery.canvas import Signature
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.annotation_versions import merge_pending_range_locks
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks import annotation_scheduler_task
from annotation.tasks.annotate_variants import annotate_variants
from annotation.tasks.annotation_scheduler_task import (
    _handle_range_lock,
    dispatch_annotation_runs,
    reclaim_stalled_annotation_runs,
)
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL = VariantAnnotationPipelineType.STRUCTURAL_VARIANT


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationDispatchTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        # A handful of real Variants (increasing pk) to anchor range-lock min/max FKs.
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(8)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    # ------------------------------------------------------------------ helpers
    def _make_lock(self, lo_idx, hi_idx, count, pipeline_types=(STANDARD,), external=False):
        lock = AnnotationRangeLock.objects.create(version=self.vav,
                                                  min_variant=self.variants[lo_idx],
                                                  max_variant=self.variants[hi_idx],
                                                  count=count)
        for pipeline_type in pipeline_types:
            AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type,
                                         external=external)
        return lock

    def _lease(self, annotation_run, attempt_count=1, expires_in=3600):
        annotation_run.leased_by = "worker"
        annotation_run.lease_expires = timezone.now() + timedelta(seconds=expires_in)
        annotation_run.attempt_count = attempt_count
        annotation_run.save()

    def _dispatch(self, slots, merge_noop=False):
        """ Run the dispatcher with a fixed pool size; optionally stub merge to isolate leasing.
            Returns the annotate_variants.apply_async mock for assertions. """
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=slots))
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            if merge_noop:
                stack.enter_context(mock.patch.object(annotation_scheduler_task,
                                                      "merge_pending_range_locks", return_value=0))
            dispatch_annotation_runs(self.vav.pk)
        return launch

    # ------------------------------------------------------------------ scheduler creates pending only
    def test_handle_range_lock_creates_pending_runs_and_dispatches_nothing(self):
        lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[0],
                                                  max_variant=self.variants[1], count=100)
        with mock.patch.object(annotate_variants, "apply_async") as launch:
            _handle_range_lock(lock)
        runs = AnnotationRun.objects.filter(annotation_range_lock=lock)
        self.assertEqual(runs.count(), 2)  # one per pipeline type
        for run in runs:
            self.assertEqual(run.status, AnnotationStatus.CREATED)
            self.assertIsNone(run.task_id)
            self.assertIsNone(run.lease_expires)
        launch.assert_not_called()

    # ------------------------------------------------------------------ capacity-limited dispatch
    def test_dispatch_launches_at_most_worker_slots(self):
        for i in range(5):
            self._make_lock(i, i, count=100)
        launch = self._dispatch(slots=2, merge_noop=True)
        self.assertEqual(launch.call_count, 2)
        leased = AnnotationRun.objects.filter(lease_expires__isnull=False)
        self.assertEqual(leased.count(), 2)
        for run in leased:
            self.assertEqual(run.attempt_count, 1)
        # Remaining 3 still pending, un-leased
        self.assertEqual(AnnotationRun.objects.filter(lease_expires__isnull=True).count(), 3)

    def test_dispatch_launches_lowest_pk_first(self):
        locks = [self._make_lock(i, i, count=100) for i in range(4)]
        launch = self._dispatch(slots=2, merge_noop=True)
        launched_run_ids = {c.args[0][0] for c in launch.call_args_list}
        expected = set(AnnotationRun.objects.filter(annotation_range_lock__in=locks[:2])
                       .values_list("pk", flat=True))
        self.assertEqual(launched_run_ids, expected)

    def test_no_capacity_launches_nothing(self):
        # One in-flight run already occupies the only slot
        lock = self._make_lock(0, 0, count=100)
        self._lease(lock.annotationrun_set.first())
        self._make_lock(1, 1, count=100)  # pending, but no capacity
        launch = self._dispatch(slots=1, merge_noop=True)
        launch.assert_not_called()

    # ------------------------------------------------------------------ latency path: free capacity, no merge
    def test_free_capacity_launches_unmerged(self):
        self._make_lock(0, 0, count=100)  # lone sub-BATCH_MIN lock
        with mock.patch.object(annotation_scheduler_task, "merge_pending_range_locks") as merge:
            launch = self._dispatch(slots=4)
        merge.assert_not_called()
        self.assertEqual(launch.call_count, 1)
        self.assertEqual(AnnotationRangeLock.objects.filter(version=self.vav).count(), 1)

    def test_merge_only_fires_when_backlog_exceeds_capacity(self):
        for i in range(3):
            self._make_lock(i, i, count=100)  # 3 dispatchable runs
        # capacity >= dispatchable -> no merge
        with mock.patch.object(annotation_scheduler_task, "merge_pending_range_locks") as merge:
            self._dispatch(slots=10)
        merge.assert_not_called()

    def test_merge_fires_when_backlog_exceeds_capacity(self):
        for i in range(3):
            self._make_lock(i, i, count=100)  # 3 dispatchable runs
        with mock.patch.object(annotation_scheduler_task, "merge_pending_range_locks",
                               return_value=0) as merge:
            self._dispatch(slots=1)
        merge.assert_called_once()

    # ------------------------------------------------------------------ merge mechanics
    def test_merge_combines_adjacent_pending_locks(self):
        for i in range(5):
            self._make_lock(i, i, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        merged = merge_pending_range_locks(self.vav, batch_max=1000)
        self.assertEqual(merged, 4)  # 5 locks -> 1 (4 absorbed)
        locks = AnnotationRangeLock.objects.filter(version=self.vav)
        self.assertEqual(locks.count(), 1)
        survivor = locks.first()
        self.assertEqual(survivor.count, 500)
        self.assertEqual(survivor.min_variant_id, self.variants[0].pk)
        self.assertEqual(survivor.max_variant_id, self.variants[4].pk)
        # Survivor keeps its own two runs; absorbed locks' runs cascade-deleted
        self.assertEqual(AnnotationRun.objects.filter(annotation_range_lock=survivor).count(), 2)
        self.assertEqual(AnnotationRun.objects.filter(annotation_range_lock__version=self.vav).count(), 2)

    def test_merge_respects_batch_max(self):
        for i in range(3):
            self._make_lock(i, i, count=400)
        merged = merge_pending_range_locks(self.vav, batch_max=1000)
        self.assertEqual(merged, 1)  # A+B=800<=1000, +C=1200>1000 stops
        locks = list(AnnotationRangeLock.objects.filter(version=self.vav).order_by("min_variant_id"))
        self.assertEqual(len(locks), 2)
        self.assertEqual(locks[0].count, 800)
        self.assertEqual(locks[1].count, 400)

    def test_merge_never_crosses_in_flight_lock(self):
        a = self._make_lock(0, 0, count=100)
        b = self._make_lock(1, 1, count=100)
        c = self._make_lock(2, 2, count=100)
        d = self._make_lock(3, 3, count=100)
        self._lease(c.annotationrun_set.first())  # C in-flight - a wall merge can't cross
        merged = merge_pending_range_locks(self.vav, batch_max=25000)
        self.assertEqual(merged, 1)  # only A absorbs B
        remaining = list(AnnotationRangeLock.objects.filter(version=self.vav).order_by("min_variant_id"))
        self.assertEqual(len(remaining), 3)  # AB, C, D
        self.assertEqual(remaining[0].pk, a.pk)
        self.assertEqual(remaining[0].count, 200)
        self.assertEqual(remaining[0].max_variant_id, self.variants[1].pk)
        self.assertEqual(remaining[1].pk, c.pk)
        self.assertEqual(remaining[2].pk, d.pk)

    # ------------------------------------------------------------------ dead-worker reclaim
    def test_expired_lease_is_reclaimed_to_created(self):
        lock = self._make_lock(0, 0, count=100)
        run = lock.annotationrun_set.first()
        self._lease(run, attempt_count=1, expires_in=-1)  # expired
        reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.CREATED)
        self.assertIsNone(run.lease_expires)
        self.assertIsNone(run.leased_by)
        self.assertEqual(run.attempt_count, 1)  # unchanged by reclaim (bumped at lease time)
        self.assertTrue(run.is_dispatchable())

    def test_live_lease_is_not_reclaimed(self):
        lock = self._make_lock(0, 0, count=100)
        run = lock.annotationrun_set.first()
        self._lease(run, attempt_count=1, expires_in=3600)  # still live
        reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertIsNotNone(run.lease_expires)
        self.assertFalse(run.is_dispatchable())

    def test_reclaim_fails_run_after_max_attempts(self):
        lock = self._make_lock(0, 0, count=100)
        run = lock.annotationrun_set.first()
        self._lease(run, attempt_count=3, expires_in=-1)  # expired, attempts exhausted (MAX=3)
        reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.ERROR)
        self.assertIsNone(run.lease_expires)
        self.assertFalse(run.is_dispatchable())

    def test_dispatch_reclaims_then_relaunches(self):
        lock = self._make_lock(0, 0, count=100)
        run = lock.annotationrun_set.first()
        self._lease(run, attempt_count=1, expires_in=-1)  # expired
        launch = self._dispatch(slots=2, merge_noop=True)
        launch.assert_called_once()
        run.refresh_from_db()
        self.assertEqual(run.attempt_count, 2)  # reclaimed -> re-leased (attempt bumped)
        self.assertIsNotNone(run.lease_expires)

    # ------------------------------------------------------------------ run completion kicks dispatcher
    def test_completion_triggers_dispatch(self):
        from annotation.tasks.annotate_variants import _trigger_dispatch
        with mock.patch.object(Signature, "apply_async") as m:
            _trigger_dispatch(self.vav.pk)
        # Immediate + 3s delayed kick (covers the commit/read race)
        self.assertEqual(m.call_count, 2)

    # ------------------------------------------------------------------ completed-but-leased run frees its slot
    def test_completed_run_does_not_count_as_in_flight(self):
        # A FINISHED run whose 4.5h lease has not been cleared must NOT shrink capacity.
        done_lock = self._make_lock(0, 0, count=100)
        done_run = done_lock.annotationrun_set.first()
        self._lease(done_run, attempt_count=1, expires_in=3600)  # live lease still set
        done_run.upload_start = timezone.now()
        done_run.upload_end = timezone.now()
        done_run.save()  # get_status() -> FINISHED
        self.assertEqual(done_run.status, AnnotationStatus.FINISHED)
        self.assertFalse(done_run.is_in_flight())

        self._make_lock(1, 1, count=100)  # one pending run
        launch = self._dispatch(slots=1, merge_noop=True)
        # The finished run's stale lease must not consume the only slot.
        launch.assert_called_once()

    # ------------------------------------------------------------------ external runs untouched
    def test_external_runs_not_dispatched(self):
        self._make_lock(0, 0, count=100, external=True)
        launch = self._dispatch(slots=4, merge_noop=True)
        launch.assert_not_called()

    def test_external_locks_not_merged(self):
        self._make_lock(0, 0, count=100, external=True)
        self._make_lock(1, 1, count=100, external=True)
        merged = merge_pending_range_locks(self.vav, batch_max=25000)
        self.assertEqual(merged, 0)
        self.assertEqual(AnnotationRangeLock.objects.filter(version=self.vav).count(), 2)

    def test_external_runs_not_reclaimed(self):
        lock = self._make_lock(0, 0, count=100, external=True)
        run = lock.annotationrun_set.first()
        self._lease(run, attempt_count=3, expires_in=-1)
        reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertNotEqual(run.status, AnnotationStatus.ERROR)
