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

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationRangeLock, AnnotationRun, AnnotationVersion, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.annotation_versions import _absorb_range_lock, merge_pending_range_locks
from annotation.tasks import annotation_scheduler_task
from annotation.tasks.annotation_scheduler_task import (
    _handle_range_lock,
    _trigger_counts_for_uncounted_runs,
    count_annotation_run,
    dispatch_annotation_runs,
    reclaim_stalled_annotation_runs,
)
from annotation.tasks.annotate_variants import annotate_variants
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
            Returns the annotate_variants.apply_async mock for assertions; the count_annotation_run
            kick mock is stored on self.count_kick. """
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=slots))
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            self.count_kick = stack.enter_context(mock.patch.object(count_annotation_run, "apply_async"))
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
        new_run = AnnotationRun.objects.create(annotation_range_lock=new_lock, pipeline_type=STANDARD)

        dispatchable_pks = {v.pk for v in
                            annotation_scheduler_task._dispatchable_variant_annotation_versions()}
        self.assertIn(new_vav.pk, dispatchable_pks)   # NEW is swept
        self.assertIn(self.vav.pk, dispatchable_pks)  # ACTIVE still swept

        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(annotation_scheduler_task, "annotation_worker_slots",
                                                  return_value=2))
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            stack.enter_context(mock.patch.object(count_annotation_run, "apply_async"))
            stack.enter_context(mock.patch.object(annotation_scheduler_task,
                                                  "merge_pending_range_locks", return_value=0))
            dispatch_annotation_runs()  # no arg -> beat sweep over all dispatchable versions

        launch.assert_any_call((new_run.pk,))
        new_run.refresh_from_db()
        self.assertEqual(new_run.attempt_count, 1)
        self.assertIsNotNone(new_run.lease_expires)

    def test_sweep_launches_uploads_before_created_across_versions(self):
        """ Global uploads-first: the no-arg sweep launches an upload-only (resume) run on ANY version
            before a full-VEP CREATED run on ANOTHER version. A CREATED run on a version iterated FIRST
            must not take the only slot ahead of an upload-only run on a version iterated later - a cheap
            DB import should never wait behind a VEP dump on another build (shared worker pool). """
        # NEW version (iterated before ACTIVE on the same build) carries only a CREATED VEP run.
        new_kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.REFSEQ, 2)
        new_kwargs["status"] = VariantAnnotationVersion.Status.NEW
        new_vav = VariantAnnotationVersion.objects.create(**new_kwargs)
        AnnotationVersion.objects.create(genome_build=self.grch37, variant_annotation_version=new_vav)
        created_lock = AnnotationRangeLock.objects.create(version=new_vav, min_variant=self.variants[0],
                                                          max_variant=self.variants[1], count=100)
        created_run = AnnotationRun.objects.create(annotation_range_lock=created_lock, pipeline_type=STANDARD)

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
                                                  return_value=1))  # only one slot - forces the choice
            launch = stack.enter_context(mock.patch.object(annotate_variants, "apply_async"))
            stack.enter_context(mock.patch.object(count_annotation_run, "apply_async"))
            stack.enter_context(mock.patch.object(annotation_scheduler_task,
                                                  "merge_pending_range_locks", return_value=0))
            dispatch_annotation_runs()  # no arg -> global sweep

        launched = {c.args[0][0] for c in launch.call_args_list}
        self.assertIn(upload_run.pk, launched)       # upload-only won the slot
        self.assertNotIn(created_run.pk, launched)   # CREATED VEP waited, despite its version coming first

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

    # ================================================================== #1646 Part 1: count + finish empties
    def test_count_task_finishes_empty_run(self):
        # The fixture holds only SNVs, so an SV run's range is empty - the count task finishes it (no dump).
        sv_run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,)).annotationrun_set.first()
        count_annotation_run(sv_run.pk)
        sv_run.refresh_from_db()
        self.assertEqual(sv_run.count, 0)
        self.assertEqual(sv_run.status, AnnotationStatus.FINISHED)
        self.assertTrue(sv_run.is_empty_finished)

    def test_count_task_fills_nonempty_count(self):
        std_run = self._make_lock(0, 2, count=100, pipeline_types=(STANDARD,)).annotationrun_set.first()
        count_annotation_run(std_run.pk)
        std_run.refresh_from_db()
        self.assertEqual(std_run.count, 3)  # variants 0,1,2 in range
        self.assertEqual(std_run.status, AnnotationStatus.CREATED)  # non-empty - stays pending

    def test_count_task_skips_leased_run(self):
        # Guarded update: a run the dispatcher already leased must not be touched by a late count task.
        sv_run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,)).annotationrun_set.first()
        self._lease(sv_run)
        count_annotation_run(sv_run.pk)
        sv_run.refresh_from_db()
        self.assertIsNone(sv_run.count)
        self.assertEqual(sv_run.status, AnnotationStatus.CREATED)

    def test_count_task_skips_run_with_no_range_lock(self):
        # #1647: the lock can be cleared after the task is kicked (reset_annotation_states, or a
        # retry-created run awaiting its lock). The count task must skip it rather than crash.
        run = self._make_lock(0, 0, count=100, pipeline_types=(STRUCTURAL,)).annotationrun_set.first()
        AnnotationRun.objects.filter(pk=run.pk).update(annotation_range_lock=None)
        count_annotation_run(run.pk)  # must not raise
        run.refresh_from_db()
        self.assertIsNone(run.count)
        self.assertEqual(run.status, AnnotationStatus.CREATED)

    def test_count_kick_targets_uncounted_runs_only(self):
        uncounted = self._make_lock(0, 0, count=100).annotationrun_set.first()
        counted = self._make_lock(1, 1, count=100).annotationrun_set.first()
        counted.count = 5
        counted.save()
        with mock.patch.object(count_annotation_run, "apply_async") as kick:
            _trigger_counts_for_uncounted_runs(self.vav)
        kicked = {c.args[0][0] for c in kick.call_args_list}
        self.assertEqual(kicked, {uncounted.pk})

    def test_dispatcher_kicks_counts_even_at_zero_capacity(self):
        # The dispatcher is the reliable count driver - it must fire even when saturated (capacity 0),
        # on db_workers, so a busy system still closes empty runs.
        uncounted = self._make_lock(0, 0, count=100).annotationrun_set.first()
        self._dispatch(slots=0, merge_noop=True)
        kicked = {c.args[0][0] for c in self.count_kick.call_args_list}
        self.assertIn(uncounted.pk, kicked)

    def test_empty_finished_run_not_dispatched(self):
        lock = self._make_lock(0, 0, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        count_annotation_run(lock.annotationrun_set.get(pipeline_type=STRUCTURAL).pk)  # SV -> empty finished
        launch = self._dispatch(slots=4, merge_noop=True)
        launch.assert_called_once()  # only the STANDARD run launches
        self.assertEqual(launch.call_args_list[0].args[0][0],
                         lock.annotationrun_set.get(pipeline_type=STANDARD).pk)

    def test_merge_combines_locks_with_empty_finished_sv_runs(self):
        # Mergeability ignores empty-finished SV runs - locks still combine on their live STANDARD runs.
        locks = [self._make_lock(i, i, count=100, pipeline_types=(STANDARD, STRUCTURAL)) for i in range(3)]
        for lock in locks:
            count_annotation_run(lock.annotationrun_set.get(pipeline_type=STRUCTURAL).pk)
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
        count_annotation_run(survivor.annotationrun_set.get(pipeline_type=STRUCTURAL).pk)  # SV -> empty finished
        absorbed = self._make_lock(1, 1, count=100, pipeline_types=(STANDARD, STRUCTURAL))
        _absorb_range_lock(survivor, absorbed)
        sv_run = survivor.annotationrun_set.get(pipeline_type=STRUCTURAL)
        std_run.refresh_from_db()
        self.assertEqual(sv_run.status, AnnotationStatus.CREATED)  # empty run reopened
        self.assertIsNone(sv_run.count)
        self.assertIsNone(std_run.count)  # stale count nulled for re-count

    # ================================================================== #1646 Part 2: resume upload-only
    def _make_past_vep_run(self, lo_idx=0, annotated_filename="/tmp/fake.vep_annotated.vcf.gz"):
        """ A run that has completed VEP (annotated VCF present) but not yet uploaded. """
        lock = self._make_lock(lo_idx, lo_idx, count=100)
        run = lock.annotationrun_set.first()
        run.dump_start = run.dump_end = timezone.now()
        run.dump_count = 100
        run.annotation_start = run.annotation_end = timezone.now()
        run.vcf_annotated_filename = annotated_filename
        run.save()
        return run

    def test_reclaim_resumes_upload_only_for_past_vep_run(self):
        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tf:
            annotated_filename = tf.name
        self.addCleanup(lambda: os.path.exists(annotated_filename) and os.remove(annotated_filename))

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

    def test_resume_run_dispatched_before_created(self):
        # A resume-upload run and a pending CREATED run, one free slot: the cheap resume run goes first.
        resume_run = self._make_past_vep_run(lo_idx=1)
        self.assertTrue(resume_run.is_upload_resumable())
        self._make_lock(0, 0, count=100)  # pending CREATED, lower pk
        launch = self._dispatch(slots=1, merge_noop=True)
        launch.assert_called_once()
        self.assertEqual(launch.call_args_list[0].args[0][0], resume_run.pk)

    def test_resume_run_not_counted_as_in_flight(self):
        # The resume run must not consume the only slot as 'in-flight' - the pending run still launches.
        self._make_past_vep_run(lo_idx=0)
        pending_run = self._make_lock(1, 1, count=100).annotationrun_set.first()
        # in-flight excludes the resume run, so capacity remains for both; with 1 slot the resume run
        # launches first. Give 2 slots so both launch, proving the resume run didn't shrink capacity.
        launch = self._dispatch(slots=2, merge_noop=True)
        self.assertEqual(launch.call_count, 2)
        launched = {c.args[0][0] for c in launch.call_args_list}
        self.assertIn(pending_run.pk, launched)

    def test_annotate_variants_resumes_upload_only(self):
        run = self._make_past_vep_run()
        self.assertEqual(run.status, AnnotationStatus.ANNOTATION_COMPLETED)
        with mock.patch("annotation.tasks.annotate_variants.dump_and_annotate_variants") as dump_mock, \
             mock.patch("annotation.tasks.annotate_variants.import_vcf_annotations") as import_mock, \
             mock.patch.object(VariantAnnotationVersion, "get_annotation_run_blocker", return_value=None), \
             mock.patch("annotation.tasks.annotate_variants.annotation_run_complete_signal"), \
             mock.patch("annotation.tasks.annotate_variants._trigger_dispatch"):
            annotate_variants.apply((run.pk,)).get()
        dump_mock.assert_not_called()  # VEP not re-run
        import_mock.assert_called_once()  # straight to upload

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
