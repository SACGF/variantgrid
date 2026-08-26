"""
Tests for the capacity-limited, merging annotation dispatcher (issue #2667).

Mirrors analysis/tests/test_scheduler.py: the scheduler now only creates pending AnnotationRuns and
dispatch_annotation_runs is the single authority that leases + launches them, merging the backlog
into bigger batches when (and only when) workers are saturated.

annotate_variants.apply_async and annotation_worker_slots are mocked so we assert dispatch decisions
without running VEP or inspecting a live celery pool.
"""
import os
import tempfile
from contextlib import ExitStack
from datetime import timedelta
from unittest import mock

from celery.canvas import Signature
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.annotation_versions import _absorb_range_lock, merge_pending_range_locks
from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    VariantAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks import annotation_scheduler_task
from annotation.annotation_run_files import get_annotated_filename
from annotation.vep_annotation import get_vep_skipped_variants_filename
from annotation.pipelines import enabled_pipeline_types, get_runner
from annotation.pipelines.vep import VEPRunner
from annotation.tasks.annotate_variants import (
    _trigger_dispatch,
    annotate_variants,
    import_annotation_run,
)
from annotation.tasks.annotation_scheduler_task import (
    _IMPORT_RUNNING_STATUSES,
    _VEP_RUNNING_STATUSES,
    COUNT_LEASE_PREFIX,
    _dispatch_counts,
    _handle_range_lock,
    _scheduled_pipeline_versions,
    _lane_in_flight_qs,
    _lease_and_launch_run,
    count_annotation_runs,
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
    def _make_lock(self, lo_idx, hi_idx, count, pipeline_types=(STANDARD,), external=False,
                   counted=True, vav=None):
        """ counted=True stamps `count` onto the runs (as the count lane would have) so the dispatcher's
            run lanes see them as ready; counted=False leaves them as count-lane candidates. """
        lock = AnnotationRangeLock.objects.create(version=vav or self.vav,
                                                  min_variant=self.variants[lo_idx],
                                                  max_variant=self.variants[hi_idx],
                                                  count=count)
        for pipeline_type in pipeline_types:
            AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=pipeline_type,
                                         external=external, count=count if counted else None)
        return lock

    def _lease(self, annotation_run, attempt_count=1, expires_in=3600):
        annotation_run.leased_by = "worker"
        annotation_run.lease_expires = timezone.now() + timedelta(seconds=expires_in)
        annotation_run.attempt_count = attempt_count
        annotation_run.save()

    def _count_run(self, annotation_run):
        """ Put a run through the count lane: null any stamped count, lease it to a count token (as
            _dispatch_counts would) and run the batched count task on it. """
        AnnotationRun.objects.filter(pk=annotation_run.pk).update(
            count=None, leased_by=f"{COUNT_LEASE_PREFIX}test",
            lease_expires=timezone.now() + timedelta(seconds=60))
        count_annotation_runs([annotation_run.pk], f"{COUNT_LEASE_PREFIX}test")

    def _dispatch(self, slots, merge_noop=False):
        """ Run the dispatcher with a fixed VEP pool size; optionally stub merge to isolate leasing.
            Returns the annotate_variants.apply_async (VEP lane) mock for assertions; the
            import_annotation_run.apply_async (import lane) mock is stored on self.import_launch and the
            count_annotation_runs kick mock on self.count_kick. """
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=slots))
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            self.import_launch = stack.enter_context(mock.patch.object(import_annotation_run, "apply_async"))
            self.count_kick = stack.enter_context(mock.patch.object(count_annotation_runs, "apply_async"))
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
            _handle_range_lock(lock, _scheduled_pipeline_versions(self.vav.genome_build))
        runs = AnnotationRun.objects.filter(annotation_range_lock=lock)
        # Only the pipelines this deployment has switched on - AnnotSV is opt-in (#720)
        self.assertEqual(runs.count(), len(enabled_pipeline_types()))
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
            self.assertEqual(run.dispatch_count, 1)
            self.assertEqual(run.attempt_count, 0)  # bumped at execution, not dispatch
        # Remaining 3 still pending, un-leased
        self.assertEqual(AnnotationRun.objects.filter(lease_expires__isnull=True).count(), 3)

    def test_no_arg_sweep_services_new_versions(self):
        """ beat: the no-arg sweep dispatches NEW versions too, not just ACTIVE. A new annotation version
            is actively built on a NEW VAV (external upload-only runs land there), so without a heartbeat
            here it only advances on discrete kicks and a stalled dispatch never self-recovers. """
        # cls.vav is ACTIVE on grch37. Add a NEW version on the same build with one pending run.
        new_kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.REFSEQ, 2)
        new_kwargs["status"] = VariantAnnotationVersion.Status.NEW
        new_vav = VariantAnnotationVersion.objects.create(**new_kwargs)
        AnnotationVersion.objects.create(genome_build=self.grch37, variant_annotation_version=new_vav)
        new_lock = AnnotationRangeLock.objects.create(version=new_vav, min_variant=self.variants[0],
                                                      max_variant=self.variants[1], count=100)
        new_run = AnnotationRun.objects.create(annotation_range_lock=new_lock, pipeline_type=STANDARD,
                                               count=100)

        dispatchable_pks = {v.pk for v in
                            annotation_scheduler_task._dispatchable_variant_annotation_versions()}
        self.assertIn(new_vav.pk, dispatchable_pks)   # NEW is swept
        self.assertIn(self.vav.pk, dispatchable_pks)  # ACTIVE still swept

        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=2))
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            stack.enter_context(mock.patch.object(count_annotation_runs, "apply_async"))
            stack.enter_context(mock.patch.object(annotation_scheduler_task,
                                                  "merge_pending_range_locks", return_value=0))
            dispatch_annotation_runs()  # no arg -> beat sweep over all dispatchable versions

        launch.assert_any_call((new_run.pk,))
        new_run.refresh_from_db()
        self.assertEqual(new_run.dispatch_count, 1)
        self.assertIsNotNone(new_run.lease_expires)

    @override_settings(ANNOTATION_UPLOAD_WORKER_SLOTS=1)
    def test_sweep_lanes_are_independent_across_versions(self):
        """ #1649: the VEP and import lanes have separate budgets in the no-arg sweep. With one VEP slot
            and one import slot, an upload-only (import lane) run on ANY version and a CREATED (VEP lane)
            run on ANOTHER version both launch in the same sweep - the cheap DB import does not consume
            the single VEP slot, and the VEP run does not consume the single import slot. """
        # NEW version (iterated before ACTIVE on the same build) carries only a CREATED VEP run.
        new_kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.REFSEQ, 2)
        new_kwargs["status"] = VariantAnnotationVersion.Status.NEW
        new_vav = VariantAnnotationVersion.objects.create(**new_kwargs)
        AnnotationVersion.objects.create(genome_build=self.grch37, variant_annotation_version=new_vav)
        created_lock = AnnotationRangeLock.objects.create(version=new_vav, min_variant=self.variants[0],
                                                          max_variant=self.variants[1], count=100)
        created_run = AnnotationRun.objects.create(annotation_range_lock=created_lock, pipeline_type=STANDARD,
                                                   count=100)

        # ACTIVE version (cls.vav, iterated after NEW) carries an upload-only run (past VEP).
        upload_lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=self.variants[2],
                                                         max_variant=self.variants[3], count=100)
        upload_run = AnnotationRun.objects.create(annotation_range_lock=upload_lock, pipeline_type=STANDARD)
        upload_run.dump_count = 100
        upload_run.annotation_start = timezone.now()
        upload_run.annotation_end = timezone.now()
        upload_run.vcf_annotated_filename = "/tmp/fake_annotated.vcf.gz"
        upload_run.save()
        self.assertEqual(upload_run.status, AnnotationStatus.ANNOTATION_COMPLETED)

        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=1))  # one VEP slot
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            import_launch = stack.enter_context(mock.patch.object(import_annotation_run, "apply_async"))
            stack.enter_context(mock.patch.object(count_annotation_runs, "apply_async"))
            stack.enter_context(mock.patch.object(annotation_scheduler_task,
                                                  "merge_pending_range_locks", return_value=0))
            dispatch_annotation_runs()  # no arg -> global sweep

        vep_launched = {c.args[0][0] for c in launch.call_args_list}
        import_launched = {c.args[0][0] for c in import_launch.call_args_list}
        self.assertEqual(import_launched, {upload_run.pk})   # import lane launched the upload-only run
        self.assertEqual(vep_launched, {created_run.pk})     # VEP lane launched the CREATED run, same sweep

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
        self.assertEqual(run.attempt_count, 1)  # re-dispatch is not an attempt (bumped at execution)
        self.assertEqual(run.dispatch_count, 1)
        self.assertIsNotNone(run.lease_expires)

    # ------------------------------------------------------------------ run completion kicks dispatcher
    def test_completion_triggers_dispatch(self):
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

    # ================================================================== #1646 Part 1: count + finish empties
    def test_count_task_finishes_empty_run(self):
        # The fixture holds only SNVs, so an SV run's range is empty - the count task finishes it (no dump).
        sv_run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,),
                                 counted=False).annotationrun_set.first()
        self._count_run(sv_run)
        sv_run.refresh_from_db()
        self.assertEqual(sv_run.count, 0)
        self.assertEqual(sv_run.status, AnnotationStatus.FINISHED)
        self.assertTrue(sv_run.is_empty_finished)
        self.assertIsNone(sv_run.leased_by)  # count lease released on the way out

    def test_count_task_fills_nonempty_count(self):
        std_run = self._make_lock(0, 2, count=100, pipeline_types=(STANDARD,),
                                  counted=False).annotationrun_set.first()
        self._count_run(std_run)
        std_run.refresh_from_db()
        self.assertEqual(std_run.count, 3)  # variants 0,1,2 in range
        self.assertEqual(std_run.status, AnnotationStatus.CREATED)  # non-empty - stays pending
        self.assertIsNone(std_run.leased_by)  # count lease released on the way out

    def test_count_task_skips_run_leased_to_someone_else(self):
        # Guarded update: a run the dispatcher leased (or another count batch holds) must not be touched
        # by a count task that no longer owns it - and its lease must be left alone.
        sv_run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,),
                                 counted=False).annotationrun_set.first()
        self._lease(sv_run)
        count_annotation_runs([sv_run.pk], f"{COUNT_LEASE_PREFIX}stale-token")
        sv_run.refresh_from_db()
        self.assertIsNone(sv_run.count)
        self.assertEqual(sv_run.status, AnnotationStatus.CREATED)
        self.assertEqual(sv_run.leased_by, "worker")  # not ours to clear

    def test_count_task_skips_run_with_no_range_lock(self):
        # #1647: the lock can be cleared after the task is kicked (reset_annotation_states, or a
        # retry-created run awaiting its lock). The count task must skip it rather than crash.
        run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,),
                              counted=False).annotationrun_set.first()
        AnnotationRun.objects.filter(pk=run.pk).update(annotation_range_lock=None)
        self._count_run(run)  # must not raise
        run.refresh_from_db()
        self.assertIsNone(run.count)
        self.assertEqual(run.status, AnnotationStatus.CREATED)
        self.assertIsNone(run.leased_by)  # lease still released

    def test_count_lane_targets_uncounted_runs_only(self):
        uncounted = self._make_lock(0, 0, count=100, counted=False).annotationrun_set.first()
        self._make_lock(1, 1, count=100)  # counted - not a count-lane candidate
        with mock.patch.object(count_annotation_runs, "apply_async") as kick:
            _dispatch_counts(self.vav, timezone.now())
        kick.assert_called_once()
        run_ids, token = kick.call_args.args[0]
        self.assertEqual(run_ids, [uncounted.pk])
        self.assertTrue(token.startswith(COUNT_LEASE_PREFIX))

    def test_count_lane_leases_batch_and_hands_out_once(self):
        # The batch is leased before it is handed out, so a second sweep in the same window (there were
        # ~250/minute in the #720 AnnotSV outage) hands out nothing - the flood is structurally impossible.
        runs = [self._make_lock(i, i, count=100, counted=False).annotationrun_set.first()
                for i in range(3)]
        now = timezone.now()
        with mock.patch.object(count_annotation_runs, "apply_async") as kick:
            _dispatch_counts(self.vav, now)
            _dispatch_counts(self.vav, now)
        kick.assert_called_once()
        run_ids, token = kick.call_args.args[0]
        self.assertEqual(set(run_ids), {run.pk for run in runs})
        for run in runs:
            run.refresh_from_db()
            self.assertEqual(run.leased_by, token)
            self.assertIsNotNone(run.lease_expires)

    def test_expired_count_lease_is_cleared_not_failed(self):
        # A count lease is not a run attempt: reclaim clears an expired one even when attempt_count is
        # exhausted, and never resets the run.
        run = self._make_lock(0, 0, count=100, counted=False).annotationrun_set.first()
        run.leased_by = f"{COUNT_LEASE_PREFIX}dead-worker"
        run.lease_expires = timezone.now() - timedelta(seconds=1)
        run.attempt_count = 3  # would be failed out if this counted as a run attempt
        run.save()
        reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.CREATED)
        self.assertIsNone(run.leased_by)
        self.assertIsNone(run.lease_expires)
        self.assertEqual(run.attempt_count, 3)  # untouched
        self.assertIsNone(run.count)  # still a count-lane candidate

    def test_count_lease_does_not_consume_vep_budget(self):
        # A count-leased batch occupies one db_worker, not VEP slots - with a whole batch out, the VEP
        # lane must still dispatch counted pending runs (the leased batch itself is held back).
        counting = self._make_lock(0, 0, count=100, counted=False).annotationrun_set.first()
        counting.leased_by = f"{COUNT_LEASE_PREFIX}in-flight"
        counting.lease_expires = timezone.now() + timedelta(seconds=60)
        counting.save()
        ready = self._make_lock(1, 1, count=100).annotationrun_set.first()
        launch = self._dispatch(slots=1, merge_noop=True)
        launch.assert_called_once()
        self.assertEqual(launch.call_args_list[0].args[0][0], ready.pk)

    def test_dispatcher_kicks_counts_even_at_zero_capacity(self):
        # The dispatcher is the reliable count driver - it must fire even when saturated (capacity 0),
        # on db_workers, so a busy system still closes empty runs.
        uncounted = self._make_lock(0, 0, count=100, counted=False).annotationrun_set.first()
        self._dispatch(slots=0, merge_noop=True)
        self.count_kick.assert_called_once()
        run_ids, _token = self.count_kick.call_args.args[0]
        self.assertIn(uncounted.pk, run_ids)

    def test_empty_finished_run_not_dispatched(self):
        lock = self._make_lock(0, 0, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        self._count_run(lock.annotationrun_set.get(pipeline_type=STRUCTURAL))  # SV -> empty finished
        launch = self._dispatch(slots=4, merge_noop=True)
        launch.assert_called_once()  # only the STANDARD run launches
        self.assertEqual(launch.call_args_list[0].args[0][0],
                         lock.annotationrun_set.get(pipeline_type=STANDARD).pk)

    def test_merge_combines_locks_with_empty_finished_sv_runs(self):
        # Mergeability ignores empty-finished SV runs - locks still combine on their live STANDARD runs.
        locks = [self._make_lock(i, i, count=100, pipeline_types=(STANDARD, STRUCTURAL)) for i in range(3)]
        for lock in locks:
            self._count_run(lock.annotationrun_set.get(pipeline_type=STRUCTURAL))
        merged = merge_pending_range_locks(self.vav, batch_max=1000)
        self.assertEqual(merged, 2)  # 3 -> 1
        survivor = AnnotationRangeLock.objects.get(version=self.vav)
        self.assertEqual(survivor.annotationrun_set.get(pipeline_type=STANDARD).status,
                         AnnotationStatus.CREATED)
        # absorb reopened the survivor's empty SV run (range grew) - now an un-counted CREATED run
        sv = survivor.annotationrun_set.get(pipeline_type=STRUCTURAL)
        self.assertEqual(sv.status, AnnotationStatus.CREATED)
        self.assertIsNone(sv.count)

    def test_merge_skips_lock_with_finished_data_run(self):
        # A lock whose STANDARD run genuinely FINISHED with data is done - not a merge participant.
        a = self._make_lock(0, 0, count=100)
        b = self._make_lock(1, 1, count=100)
        c = self._make_lock(2, 2, count=100)
        run_a = a.annotationrun_set.get(pipeline_type=STANDARD)
        run_a.dump_start = run_a.dump_end = timezone.now()
        run_a.dump_count = 100
        run_a.annotation_start = run_a.annotation_end = timezone.now()
        run_a.upload_start = run_a.upload_end = timezone.now()
        run_a.save()  # FINISHED with data (dump_count > 0) -> not is_empty_finished
        self.assertFalse(run_a.is_empty_finished)
        merged = merge_pending_range_locks(self.vav, batch_max=25000)
        self.assertEqual(merged, 1)  # b absorbs c; a left alone
        remaining = list(AnnotationRangeLock.objects.filter(version=self.vav).order_by("min_variant_id"))
        self.assertEqual([lock.pk for lock in remaining], [a.pk, b.pk])
        self.assertEqual(remaining[1].max_variant_id, self.variants[2].pk)

    def test_absorb_reopens_empty_finished_run_and_nulls_counts(self):
        survivor = self._make_lock(0, 0, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        std_run = survivor.annotationrun_set.get(pipeline_type=STANDARD)
        std_run.count = 1  # counted non-empty
        std_run.save()
        self._count_run(survivor.annotationrun_set.get(pipeline_type=STRUCTURAL))  # SV -> empty finished
        absorbed = self._make_lock(1, 1, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        _absorb_range_lock(survivor, absorbed)
        sv_run = survivor.annotationrun_set.get(pipeline_type=STRUCTURAL)
        std_run.refresh_from_db()
        self.assertEqual(sv_run.status, AnnotationStatus.CREATED)  # empty run reopened
        self.assertIsNone(sv_run.count)
        self.assertIsNone(std_run.count)  # stale count nulled for re-count

    # ================================================================== #1646 Part 2: resume upload-only
    def _make_past_vep_run(self, lo_idx=0, annotated_filename="/tmp/fake.vep_annotated.vcf.gz", vav=None):
        """ A run that has completed VEP (annotated VCF present) but not yet uploaded. """
        lock = self._make_lock(lo_idx, lo_idx, count=100, vav=vav)
        run = lock.annotationrun_set.first()
        run.dump_start = run.dump_end = timezone.now()
        run.dump_count = 100
        run.annotation_start = run.annotation_end = timezone.now()
        run.vcf_annotated_filename = annotated_filename
        run.save()
        return run

    def _real_annotated_file(self):
        """ An annotated VCF that actually exists - import_annotation_run re-checks the file on disk at
            execution time, so tests that execute the import task need a real one. """
        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tf:
            annotated_filename = tf.name
        self.addCleanup(lambda: os.path.exists(annotated_filename) and os.remove(annotated_filename))
        return annotated_filename

    def test_reclaim_resumes_upload_only_for_past_vep_run(self):
        annotated_filename = self._real_annotated_file()
        run = self._make_past_vep_run(annotated_filename=annotated_filename)
        run.upload_start = timezone.now()  # died mid-upload -> UPLOAD_STARTED
        run.save()
        self._lease(run, attempt_count=1, expires_in=-1)  # expired lease
        self.assertEqual(run.status, AnnotationStatus.UPLOAD_STARTED)

        with mock.patch.object(AnnotationRun, "delete_related_objects") as dro:
            reclaim_stalled_annotation_runs(self.vav)
        dro.assert_called_once()  # partial upload rows scrubbed
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        self.assertEqual(run.vcf_annotated_filename, annotated_filename)
        self.assertTrue(os.path.exists(annotated_filename))  # VEP output kept
        self.assertIsNone(run.upload_start)
        self.assertIsNone(run.lease_expires)
        self.assertIsNone(run.task_id)
        self.assertTrue(run.is_upload_resumable())
        self.assertFalse(run.is_in_flight())

    def test_resume_run_dispatched_on_import_lane(self):
        # #1649: a resume/upload-only run (ANNOTATION_COMPLETED, past VEP) launches via import_annotation_run
        # (db_workers), while a pending CREATED run launches via annotate_variants (VEP). Independent lanes.
        resume_run = self._make_past_vep_run(lo_idx=1)
        self.assertTrue(resume_run.is_upload_resumable())
        created_run = self._make_lock(0, 0, count=100).annotationrun_set.first()  # pending CREATED, lower pk
        launch = self._dispatch(slots=1, merge_noop=True)
        self.import_launch.assert_called_once()
        self.assertEqual(self.import_launch.call_args_list[0].args[0][0], resume_run.pk)
        launch.assert_called_once()
        self.assertEqual(launch.call_args_list[0].args[0][0], created_run.pk)

    def test_resume_run_does_not_shrink_vep_capacity(self):
        # #1649: an ANNOTATION_COMPLETED resume run lives in the import lane, so it must not count against
        # the VEP budget. With one VEP slot the pending CREATED run still launches on the VEP lane, and the
        # resume run launches on the import lane in the same cycle.
        resume_run = self._make_past_vep_run(lo_idx=0)
        pending_run = self._make_lock(1, 1, count=100).annotationrun_set.first()
        launch = self._dispatch(slots=1, merge_noop=True)
        launch.assert_called_once()
        self.assertEqual(launch.call_args_list[0].args[0][0], pending_run.pk)
        self.import_launch.assert_called_once()
        self.assertEqual(self.import_launch.call_args_list[0].args[0][0], resume_run.pk)

    def test_annotate_variants_is_vep_only(self):
        # #1649: annotate_variants dumps + runs VEP and stops at ANNOTATION_COMPLETED. It does NOT import
        # (the dispatcher runs import_annotation_run on db_workers) and does NOT send the complete signal.
        run = self._make_lock(0, 0, count=100).annotationrun_set.first()
        self.assertEqual(run.status, AnnotationStatus.CREATED)
        with mock.patch.object(VEPRunner, "annotate") as dump_mock, \
             mock.patch.object(VEPRunner, "import_results") as import_mock, \
             mock.patch.object(VariantAnnotationVersion, "get_annotation_run_blocker", return_value=None), \
             mock.patch("annotation.tasks.annotate_variants.annotation_run_complete_signal") as signal, \
             mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
            annotate_variants.apply((run.pk,)).get()
        dump_mock.assert_called_once()   # VEP ran
        import_mock.assert_not_called()  # upload deferred to the import lane
        signal.send.assert_not_called()  # not complete until imported

    def test_import_annotation_run_imports_and_signals(self):
        # #1649: import_annotation_run bulk-loads the annotated VCF and sends the complete signal once.
        run = self._make_past_vep_run(annotated_filename=self._real_annotated_file())
        self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        with mock.patch.object(VEPRunner, "import_results") as import_mock, \
             mock.patch("annotation.tasks.annotate_variants.annotation_run_complete_signal") as signal, \
             mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
            import_annotation_run.apply((run.pk,)).get()
        import_mock.assert_called_once()   # straight to upload
        signal.send.assert_called_once()   # complete signal fires from the import lane
        run.refresh_from_db()
        self.assertIsNone(run.task_id)     # lease/lock released in finally
        self.assertIsNone(run.lease_expires)

    def test_reclaim_full_scrub_when_no_annotated_file(self):
        # Pre-VEP progress (no annotated VCF): existing behaviour - full scrub back to CREATED.
        lock = self._make_lock(0, 0, count=100)
        run = lock.annotationrun_set.first()
        run.dump_start = timezone.now()
        run.dump_count = 100
        run.vcf_dump_filename = None
        run.save()
        self._lease(run, attempt_count=1, expires_in=-1)
        with mock.patch.object(AnnotationRun, "delete_related_objects"):
            reclaim_stalled_annotation_runs(self.vav)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.CREATED)
        self.assertIsNone(run.dump_start)
        self.assertIsNone(run.dump_count)
        self.assertTrue(run.is_dispatchable())

    # ================================================================== #1660: reclaim derives paths
    def _run_with_dump(self, tmp_dir, lo_idx=0):
        """ A run leased, expired, and part-way through the pipeline with its dump stem persisted -
            which is all the row carries until dump_and_annotate_variants' final save. """
        run = self._make_lock(lo_idx, lo_idx, count=100).annotationrun_set.first()
        run.dump_start = timezone.now()
        run.dump_end = timezone.now()
        run.dump_count = 100
        run.annotation_start = timezone.now()
        run.vcf_dump_filename = os.path.join(tmp_dir, f"dump_{run.pk}__task.vcf.gz")
        run.save()
        self._lease(run, attempt_count=1, expires_in=-1)
        return run

    def test_reclaim_resumes_upload_only_when_row_has_no_annotated_filename(self):
        # #1660: vcf_annotated_filename only reaches the DB at the final save, after VEP *and* AnnotSV.
        # A run reclaimed mid-AnnotSV has a complete annotated VCF on disk but a NULL field, so reclaim
        # must derive the path from the dump stem - reading the field would scrub finished VEP work.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run = self._run_with_dump(tmp_dir)
            annotated_filename = get_annotated_filename(run, run.vcf_dump_filename)
            skipped_variants_filename = get_vep_skipped_variants_filename(annotated_filename)
            open(annotated_filename, "w").close()  # VEP finished; row not saved yet
            open(skipped_variants_filename, "w").close()  # ... written by Runner::finish()
            self.assertIsNone(run.vcf_annotated_filename)

            with mock.patch.object(AnnotationRun, "delete_related_objects"):
                reclaim_stalled_annotation_runs(self.vav)

            run.refresh_from_db()
            self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
            self.assertTrue(os.path.exists(annotated_filename))  # VEP output kept
            # Derived paths recorded, so the upload-only relaunch can find what we kept
            self.assertEqual(run.vcf_annotated_filename, annotated_filename)
            self.assertEqual(run.vep_skipped_variants_filename, skipped_variants_filename)
            self.assertTrue(run.is_upload_resumable())

    def test_reclaim_scrubs_run_killed_mid_vep(self):
        # #1710: VEP creates the annotated VCF at startup and writes to it progressively, so a run killed
        # mid-VEP leaves a part-written file. Resuming that upload-only sent a truncated VCF to the import
        # lane, where the #1701 checks failed it - turning one stalled run into an ERROR. Only the
        # skipped-variants list (written by Runner::finish, after the output handle is closed) says VEP
        # got to the end, so without it the run goes back for a full re-dump.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run = self._run_with_dump(tmp_dir)
            annotated_filename = get_annotated_filename(run, run.vcf_dump_filename)
            with open(annotated_filename, "w") as f:
                f.write("##fileformat=VCFv4.2\n")  # VEP got as far as the header
            self.assertFalse(os.path.exists(get_vep_skipped_variants_filename(annotated_filename)))

            with mock.patch.object(AnnotationRun, "delete_related_objects"):
                reclaim_stalled_annotation_runs(self.vav)

            run.refresh_from_db()
            self.assertEqual(run.status, AnnotationStatus.CREATED)
            self.assertTrue(run.is_dispatchable())
            self.assertFalse(run.is_upload_resumable())
            self.assertIsNone(run.vcf_annotated_filename)
            self.assertIsNone(run.annotation_end)
            self.assertIsNone(run.dump_count)  # re-dumped from scratch, so the old count can't stand
            # The part-written VCF is collected, so the re-dump doesn't land next to it
            self.assertFalse(os.path.exists(annotated_filename))

    def test_reclaim_full_scrub_removes_derived_files(self):
        # #1660: the row names only the dump; everything else the runner writes is derivable from that
        # stem. Reclaim must collect them - a genuinely dead worker never runs
        # _cleanup_reclaimed_run_files, so reclaim is the only collector for that case.
        with tempfile.TemporaryDirectory() as tmp_dir, \
                override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
            run = self._run_with_dump(tmp_dir)
            dump_filename = run.vcf_dump_filename
            annotated_filename = get_annotated_filename(run, dump_filename)
            # Partial VEP output - the heartbeat's kill leaves this behind
            derived = get_runner(run.pipeline_type).get_output_paths(run, dump_filename)
            self.assertIn(annotated_filename, derived)
            for path in derived:
                open(path, "w").close()
            # Reclaim takes the scrub branch only when there's no annotated VCF to resume from
            os.remove(annotated_filename)

            with mock.patch.object(AnnotationRun, "delete_related_objects"):
                reclaim_stalled_annotation_runs(self.vav)

            for path in derived:
                self.assertFalse(os.path.exists(path), f"reclaim left {os.path.basename(path)} behind")
            run.refresh_from_db()
            self.assertEqual(run.status, AnnotationStatus.CREATED)
            self.assertIsNone(run.vcf_dump_filename)
            self.assertTrue(run.is_dispatchable())

    # ================================================================== #1649: VEP / import lane split
    def test_dispatch_routes_by_status(self):
        # CREATED -> annotate_variants (VEP lane); ANNOTATION_COMPLETED-with-file -> import_annotation_run.
        created_run = self._make_lock(0, 0, count=100).annotationrun_set.first()
        import_run = self._make_past_vep_run(lo_idx=1)
        launch = self._dispatch(slots=4, merge_noop=True)
        vep_launched = {c.args[0][0] for c in launch.call_args_list}
        import_launched = {c.args[0][0] for c in self.import_launch.call_args_list}
        self.assertEqual(vep_launched, {created_run.pk})
        self.assertEqual(import_launched, {import_run.pk})

    @override_settings(ANNOTATION_UPLOAD_WORKER_SLOTS=4)
    def test_saturated_vep_lane_does_not_block_import(self):
        # A full VEP lane (all VEP slots leased) must not stop the import lane launching an upload-only run.
        for i in range(2):  # fill both VEP slots with in-flight CREATED runs
            self._lease(self._make_lock(i, i, count=100).annotationrun_set.first())
        import_run = self._make_past_vep_run(lo_idx=3)
        # Another CREATED run that can't launch (VEP saturated)
        blocked_created = self._make_lock(4, 4, count=100).annotationrun_set.first()
        launch = self._dispatch(slots=2, merge_noop=True)
        launch.assert_not_called()  # VEP lane saturated - nothing new on it
        self.import_launch.assert_called_once()
        self.assertEqual(self.import_launch.call_args_list[0].args[0][0], import_run.pk)
        blocked_created.refresh_from_db()
        self.assertIsNone(blocked_created.lease_expires)  # stayed pending

    @override_settings(ANNOTATION_UPLOAD_WORKER_SLOTS=1)
    def test_saturated_import_lane_does_not_block_vep(self):
        # A full import lane (its one slot leased) must not stop the VEP lane launching a CREATED run.
        busy_import = self._make_past_vep_run(lo_idx=0)
        self._lease(busy_import)  # occupies the single import slot
        created_run = self._make_lock(1, 1, count=100).annotationrun_set.first()
        another_import = self._make_past_vep_run(lo_idx=2)  # can't launch - import lane full
        launch = self._dispatch(slots=2, merge_noop=True)
        self.import_launch.assert_not_called()  # import lane saturated
        launch.assert_called_once()
        self.assertEqual(launch.call_args_list[0].args[0][0], created_run.pk)
        another_import.refresh_from_db()
        self.assertIsNone(another_import.lease_expires)  # stayed pending

    def test_lane_in_flight_accounting(self):
        # A leased VEP-phase run counts only in the VEP lane; a leased import-phase run only in the import
        # lane. Neither leaks into the other lane's in-flight budget.
        vep_run = self._make_lock(0, 0, count=100).annotationrun_set.first()  # CREATED
        self._lease(vep_run)
        import_run = self._make_past_vep_run(lo_idx=1)  # ANNOTATION_COMPLETED
        self._lease(import_run)
        now = timezone.now()
        vep_in_flight = set(_lane_in_flight_qs(self.vav, now, _VEP_RUNNING_STATUSES)
                            .values_list("pk", flat=True))
        import_in_flight = set(_lane_in_flight_qs(self.vav, now, _IMPORT_RUNNING_STATUSES)
                               .values_list("pk", flat=True))
        self.assertEqual(vep_in_flight, {vep_run.pk})
        self.assertEqual(import_in_flight, {import_run.pk})

    def test_pending_runs_not_counted_in_flight_in_either_lane(self):
        # Un-leased CREATED (VEP lane) and un-leased ANNOTATION_COMPLETED (import lane) hold no slot yet.
        self._make_lock(0, 0, count=100)                # pending CREATED
        self._make_past_vep_run(lo_idx=1)               # un-leased ANNOTATION_COMPLETED
        now = timezone.now()
        self.assertEqual(_lane_in_flight_qs(self.vav, now, _VEP_RUNNING_STATUSES).count(), 0)
        self.assertEqual(_lane_in_flight_qs(self.vav, now, _IMPORT_RUNNING_STATUSES).count(), 0)

    def test_external_imported_run_launches_on_import_lane(self):
        # #1568/#1649: an external run whose external VEP was imported (external cleared, ANNOTATION_COMPLETED
        # with an annotated VCF) rejoins post-VEP and launches on the import lane, never annotate_variants.
        lock = self._make_lock(0, 0, count=100, external=True)
        run = lock.annotationrun_set.first()
        # Mirror _import_annotated_annotation_run: clear external, stamp annotation dates + annotated file.
        run.external = False
        run.dump_count = 100
        run.annotation_start = run.annotation_end = timezone.now()
        run.vcf_annotated_filename = "/tmp/fake_external.vep_annotated.vcf.gz"
        run.save()
        self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        launch = self._dispatch(slots=4, merge_noop=True)
        launch.assert_not_called()  # never the VEP lane
        self.import_launch.assert_called_once()
        self.assertEqual(self.import_launch.call_args_list[0].args[0][0], run.pk)

    def test_reclaimed_import_run_relaunches_on_import_lane(self):
        # A stalled import-phase run is reclaimed to ANNOTATION_COMPLETED (annotated file kept) and
        # re-launched on the import lane next cycle - not the VEP lane.
        run = self._make_past_vep_run(annotated_filename=self._real_annotated_file())
        run.upload_start = timezone.now()  # died mid-import -> UPLOAD_STARTED
        run.save()
        self._lease(run, attempt_count=1, expires_in=-1)  # expired lease
        with mock.patch.object(AnnotationRun, "delete_related_objects"):
            launch = self._dispatch(slots=4, merge_noop=True)
        run.refresh_from_db()
        self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        launch.assert_not_called()  # not the VEP lane
        self.import_launch.assert_called_once()
        self.assertEqual(self.import_launch.call_args_list[0].args[0][0], run.pk)
        self.assertEqual(run.attempt_count, 1)  # re-dispatch is not an attempt (bumped at execution)

    # ================================================================== over-committed dispatch fixes
    def _make_new_vav(self):
        """ A second dispatchable version on the same build (latest NEW, alongside the ACTIVE cls.vav). """
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.REFSEQ, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.NEW
        new_vav = VariantAnnotationVersion.objects.create(**kwargs)
        AnnotationVersion.objects.create(genome_build=self.grch37, variant_annotation_version=new_vav)
        return new_vav

    @override_settings(ANNOTATION_UPLOAD_WORKER_SLOTS=2)
    def test_per_vav_kick_budgets_import_capacity_globally(self):
        # The bug behind the AnnotSV import-lane stall: the per-completion kick budgeted capacity per-VAV,
        # so each version saw the full pool minus only its own in-flight and collectively over-committed
        # it. With the global budget already filled by another version's in-flight imports, a kick for
        # this version must launch nothing.
        other_vav = self._make_new_vav()
        for i in range(2):  # fill the whole global import budget from the other version
            self._lease(self._make_past_vep_run(lo_idx=i, vav=other_vav))
        self._make_past_vep_run(lo_idx=3)  # waiting on cls.vav
        self._dispatch(slots=4, merge_noop=True)  # per-completion path: dispatch_annotation_runs(vav.pk)
        self.import_launch.assert_not_called()

    def test_lease_and_launch_run_is_atomic(self):
        # The conditional lease loses (returns False, launches nothing) when the run was leased, locked,
        # or changed status between the dispatcher's read and the lease write.
        now = timezone.now()
        with mock.patch.object(annotate_variants, "apply_async") as launch, \
                mock.patch.object(import_annotation_run, "apply_async") as import_launch:
            already_leased = self._make_lock(0, 0, count=100).annotationrun_set.first()
            self._lease(already_leased)
            self.assertFalse(_lease_and_launch_run(already_leased, "dispatch:test", now))

            locked = self._make_lock(1, 1, count=100).annotationrun_set.first()
            AnnotationRun.objects.filter(pk=locked.pk).update(task_id="some-task")
            self.assertFalse(_lease_and_launch_run(locked, "dispatch:test", now))

            moved_on = self._make_lock(2, 2, count=100).annotationrun_set.first()  # read as CREATED...
            AnnotationRun.objects.filter(pk=moved_on.pk).update(status=AnnotationStatus.FINISHED)
            self.assertFalse(_lease_and_launch_run(moved_on, "dispatch:test", now))

            launch.assert_not_called()
            import_launch.assert_not_called()

            fresh = self._make_lock(3, 3, count=100).annotationrun_set.first()
            self.assertTrue(_lease_and_launch_run(fresh, "dispatch:test", now))
            launch.assert_called_once()
        fresh.refresh_from_db()
        self.assertEqual(fresh.dispatch_count, 1)
        self.assertEqual(fresh.attempt_count, 0)

    def test_stale_import_task_is_noop_on_finished_run(self):
        # Shape 1 of the AnnotSV outage: the run was double-dispatched, the first import succeeded and
        # cleanup removed the annotated file. The second, still-queued task must be a silent no-op -
        # no FileNotFoundError, no error_exception - with the lock released.
        run = self._make_past_vep_run(annotated_filename="/does/not/exist.vep_annotated.vcf.gz")
        run.upload_start = run.upload_end = timezone.now()
        run.annotated_count = 100
        run.save()
        self.assertEqual(run.status, AnnotationStatus.FINISHED)
        with mock.patch.object(VEPRunner, "import_results") as import_mock, \
             mock.patch("annotation.tasks.annotate_variants.annotation_run_complete_signal") as signal, \
             mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
            import_annotation_run.apply((run.pk,)).get()
        import_mock.assert_not_called()   # no second import
        signal.send.assert_not_called()
        run.refresh_from_db()
        self.assertIsNone(run.error_exception)
        self.assertEqual(run.status, AnnotationStatus.FINISHED)
        self.assertIsNone(run.task_id)    # lock released
        self.assertIsNone(run.lease_expires)

    def test_attempt_count_bumped_by_execution_not_dispatch(self):
        run = self._make_lock(0, 0, count=100).annotationrun_set.first()
        launch = self._dispatch(slots=1, merge_noop=True)
        launch.assert_called_once()
        run.refresh_from_db()
        self.assertEqual(run.dispatch_count, 1)
        self.assertEqual(run.attempt_count, 0)  # queued, not yet executed - no retry budget burned
        with mock.patch.object(VEPRunner, "annotate"), \
             mock.patch.object(VariantAnnotationVersion, "get_annotation_run_blocker", return_value=None), \
             mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
            annotate_variants.apply((run.pk,)).get()
        run.refresh_from_db()
        self.assertEqual(run.attempt_count, 1)  # bumped exactly once, when the task took its lock
