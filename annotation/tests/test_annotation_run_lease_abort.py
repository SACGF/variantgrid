"""
Tests for the reclaimed-run safety net (issue #1658): when a stalled worker loses its lease and the
dispatcher hands the AnnotationRun to a fresh attempt, the losing attempt must

  * write its VEP output to a path the winner is not also writing to (unique per-task dump paths), and
  * abort promptly (the lease heartbeat kills the VEP subprocess), and
  * unwind without touching a row it no longer owns.

The first two are covered by e707e13f5. The third is not - see AnnotationRunUnwindOwnershipTest.

These deliberately avoid asserting the *format* of the dump path: the invariant that matters is
"two concurrent attempts never collide, and external dumps (#1568) stay stable", not the exact stem.
"""
import os
import tempfile
import threading
import time
from datetime import timedelta
from unittest import mock

from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import AnnotationRangeLock, AnnotationRun, AnnotationVersion, VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.tasks import annotate_variants as annotate_variants_module
from annotation.tasks.annotate_variants import (
    AnnotationRunLeaseHeartbeat,
    annotate_variants,
    dump_and_annotate_variants,
    import_annotation_run,
)
from genes.models_enums import AnnotationConsortium
from library.utils import execute_cmd
from library.utils.file_utils import name_from_filename
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

STANDARD = VariantAnnotationPipelineType.STANDARD


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class AnnotationRunLeaseTestCase(TestCase):
    """ Shared fixture: one ACTIVE VariantAnnotationVersion and a helper to make leased AnnotationRuns. """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.variants = [slowly_create_test_variant("1", 100000 + i * 10, 'A', 'T', cls.grch37)
                        for i in range(2)]
        kwargs = get_fake_vep_version(cls.grch37, AnnotationConsortium.ENSEMBL, 2)
        kwargs["status"] = VariantAnnotationVersion.Status.ACTIVE
        cls.vav = VariantAnnotationVersion.objects.create(**kwargs)
        cls.av = AnnotationVersion.objects.create(genome_build=cls.grch37, variant_annotation_version=cls.vav)

    def setUp(self):
        super().setUp()
        self._dump_dir = tempfile.mkdtemp(prefix="annotation_dump_")
        patcher = override_settings(ANNOTATION_VCF_DUMP_DIR=self._dump_dir)
        patcher.enable()
        self.addCleanup(patcher.disable)

    def _make_run(self, task_id=None, external=False) -> AnnotationRun:
        lock = AnnotationRangeLock.objects.create(version=self.vav,
                                                  min_variant=self.variants[0],
                                                  max_variant=self.variants[-1],
                                                  count=10)
        run = AnnotationRun.objects.create(annotation_range_lock=lock, pipeline_type=STANDARD,
                                           external=external, count=10)
        if task_id:
            run.task_id = task_id
            run.leased_by = "worker-a"
            run.lease_expires = timezone.now() + timedelta(seconds=600)
            run.save()
        return run


class AnnotationRunUnwindOwnershipTest(AnnotationRunLeaseTestCase):
    """ #1658: the losing attempt's try/except/finally in annotate_variants operates on an ORM instance
        loaded before the reclaim, and saves it whole. Once the run has been handed to a new attempt,
        every one of those saves is writing to a row this worker no longer owns. """

    def _run_losing_attempt(self, reclaim_side_effect):
        """ Execute annotate_variants eagerly (so it gets a real task_id) with the VEP stage replaced by
            `reclaim_side_effect`, which simulates the run being reclaimed mid-VEP and this attempt then
            being killed by its own lease heartbeat. Returns the EagerResult. """
        run = self._make_run()
        with mock.patch.object(annotate_variants_module, "dump_and_annotate_variants",
                               side_effect=reclaim_side_effect), \
                mock.patch.object(annotate_variants_module, "_trigger_dispatch"), \
                mock.patch.object(annotate_variants_module, "report_message"), \
                mock.patch.object(annotate_variants_module, "create_event"), \
                mock.patch.object(VariantAnnotationVersion, "get_annotation_run_blocker", return_value=None):
            result = annotate_variants.apply(args=[run.pk])
        return run, result

    def test_losing_attempt_unwind_leaves_new_owners_lease_intact(self):
        """ The `finally` in annotate_variants clears task_id/leased_by/lease_expires and saves. If the run
            was reclaimed while we were in VEP, those fields belong to the *new* attempt - blanking them
            re-opens the run to the dispatcher while the new attempt is still running VEP, and the new
            attempt's own heartbeat then sees a task_id mismatch and aborts itself. One reclaim cascades. """
        new_lease_expires = timezone.now() + timedelta(seconds=600)

        def _reclaimed_then_killed(annotation_run, **kwargs):
            AnnotationRun.objects.filter(pk=annotation_run.pk).update(
                task_id="new-owner-task", leased_by="worker-b", lease_expires=new_lease_expires)
            raise RuntimeError("VEP returned -9")  # heartbeat killed our subprocess

        run, result = self._run_losing_attempt(_reclaimed_then_killed)
        self.assertTrue(result.failed())

        run.refresh_from_db()
        self.assertEqual(run.task_id, "new-owner-task",
                         "losing attempt cleared the new owner's execution lock")
        self.assertEqual(run.leased_by, "worker-b",
                         "losing attempt cleared the new owner's lease holder")
        self.assertIsNotNone(run.lease_expires,
                             "losing attempt released the new owner's lease - run is now re-dispatchable "
                             "while the new owner is still running VEP")

    def test_losing_attempt_unwind_does_not_write_stale_fields(self):
        """ The same stale-instance saves also push this attempt's own dump path and traceback onto the
            new owner's row: vcf_dump_filename then points at the loser's file (which the loser may have
            left incomplete) while the winner writes somewhere else, and error_exception flips the winner
            to ERROR status mid-run. """
        winner_dump = os.path.join(self._dump_dir, "winner__task-new-owner.vcf")

        def _reclaimed_then_killed(annotation_run, **kwargs):
            # This attempt got as far as dumping before being reclaimed.
            annotation_run.vcf_dump_filename = os.path.join(self._dump_dir, "loser__task-old.vcf")
            AnnotationRun.objects.filter(pk=annotation_run.pk).update(
                task_id="new-owner-task", leased_by="worker-b",
                lease_expires=timezone.now() + timedelta(seconds=600),
                vcf_dump_filename=winner_dump)
            raise RuntimeError("VEP returned -9")

        run, result = self._run_losing_attempt(_reclaimed_then_killed)
        self.assertTrue(result.failed())

        run.refresh_from_db()
        self.assertEqual(run.vcf_dump_filename, winner_dump,
                         "losing attempt overwrote the new owner's dump path with its own")
        self.assertIsNone(run.error_exception,
                          "losing attempt's traceback was written onto the new owner's run")

    def test_losing_attempt_removes_its_own_orphaned_output_files(self):
        """ Per-task dump paths (#1658) mean the losing attempt's files are named after *its* task_id. The
            reclaim deletes only what the DB row names and then NULLs those fields, so once the run is
            handed on, nothing in the pipeline can ever find this attempt's dump/annotated VCF again -
            there is no sweep of ANNOTATION_VCF_DUMP_DIR. The losing attempt is the last party that knows
            the paths, so it must delete them itself or every reclaim leaks a VCF pair to disk. """
        loser_dump = os.path.join(self._dump_dir, "loser__task-old.vcf")
        # Names the annotated VCF exactly as dump_and_annotate_variants derives it from the dump.
        loser_annotated = os.path.join(
            self._dump_dir, f"{name_from_filename(loser_dump)}.vep_annotated_{self.grch37.name}.vcf.gz")
        for path in (loser_dump, loser_annotated):
            with open(path, "wt") as f:
                f.write("partial VEP output\n")

        def _reclaimed_then_killed(annotation_run, **kwargs):
            annotation_run.vcf_dump_filename = loser_dump
            # Reclaim: deletes/NULLs what the row names (the row named nothing yet) and re-leases.
            AnnotationRun.objects.filter(pk=annotation_run.pk).update(
                task_id="new-owner-task", leased_by="worker-b", vcf_dump_filename=None,
                lease_expires=timezone.now() + timedelta(seconds=600))
            raise RuntimeError("VEP returned -9")

        run, result = self._run_losing_attempt(_reclaimed_then_killed)
        self.assertTrue(result.failed())

        self.assertFalse(os.path.exists(loser_dump),
                         "losing attempt left its dump VCF orphaned on disk - unreachable by any later "
                         "reclaim, which only deletes what the DB row names")
        self.assertFalse(os.path.exists(loser_annotated),
                         "losing attempt left its (possibly partial) annotated VCF orphaned on disk")


class AnnotationRunLockFailureTest(AnnotationRunLeaseTestCase):
    """ #1658: the task_id lock-failure branch is the one place in the pipeline that has *proved* it does
        not own the run - it is only reached because another task won the lock. It must therefore write
        nothing but its own audit entry. This is the same clobber as the reclaim path, minus the reclaim:
        an ordinary double-dispatch is enough to trip it. """

    def _stale_instance_loaded_before_the_winner(self, run) -> AnnotationRun:
        """ Model the TOCTOU window the branch exists for: this task loaded the run while it was free and
            the winner claimed it before our conditional UPDATE ran, so the instance we hold still says
            task_id=None - exactly the value a full save() would put over the winner's lock. """
        stale = AnnotationRun.objects.get(pk=run.pk)
        stale.task_id = None
        stale.leased_by = None
        stale.lease_expires = None
        return stale

    def _assert_winners_lock_intact(self, run, lease_expires):
        run.refresh_from_db()
        self.assertEqual(run.task_id, "winner-task",
                         "lock loser cleared the winning task's execution lock")
        self.assertEqual(run.leased_by, "worker-a",
                         "lock loser cleared the winning task's lease holder")
        self.assertEqual(run.lease_expires, lease_expires,
                         "lock loser released the winning task's lease - the run is now re-dispatchable "
                         "while the winner is still running")

    def test_annotate_lock_failure_leaves_the_winners_lock_intact(self):
        run = self._make_run(task_id="winner-task")
        lease_expires = run.lease_expires
        stale = self._stale_instance_loaded_before_the_winner(run)

        with mock.patch.object(AnnotationRun.objects, "get", return_value=stale):
            result = annotate_variants.apply(args=[run.pk])
        self.assertTrue(result.failed())

        self._assert_winners_lock_intact(run, lease_expires)
        self.assertIn("lock_failed", run.celery_task_logs.get(result.id, {}),
                      "lock failure was not recorded against the losing task")
        self.assertNotIn("lock_failed", run.celery_task_logs.get("winner-task", {}),
                         "lock failure was recorded against the winning task")

    def test_import_lock_failure_leaves_the_winners_lock_intact(self):
        """ import_annotation_run carries a copy of the same lock-acquisition block, so it needs the same
            guard - and its own test, or the two drift. """
        run = self._make_run(task_id="winner-task")
        lease_expires = run.lease_expires
        stale = self._stale_instance_loaded_before_the_winner(run)

        with mock.patch.object(AnnotationRun.objects, "get", return_value=stale):
            result = import_annotation_run.apply(args=[run.pk])
        self.assertTrue(result.failed())

        self._assert_winners_lock_intact(run, lease_expires)
        self.assertIn("lock_failed", run.celery_task_logs.get(result.id, {}),
                      "lock failure was not recorded against the losing task")

    def test_lock_failure_keeps_the_winners_concurrent_task_log(self):
        """ Recording the failure is a read-modify-write of celery_task_logs. Doing it from the instance
            we loaded before the winner existed would silently drop whatever the winner has logged since. """
        run = self._make_run(task_id="winner-task")
        stale = self._stale_instance_loaded_before_the_winner(run)
        # Winner logs its start after our instance was loaded.
        run.celery_task_logs = {"winner-task": {"start": "2026-07-21T00:00:00"}}
        run.save()

        with mock.patch.object(AnnotationRun.objects, "get", return_value=stale):
            result = annotate_variants.apply(args=[run.pk])
        self.assertTrue(result.failed())

        run.refresh_from_db()
        self.assertEqual(run.celery_task_logs.get("winner-task", {}).get("start"), "2026-07-21T00:00:00",
                         "lock loser overwrote celery_task_logs from its stale copy, losing the winner's entry")
        self.assertIn("lock_failed", run.celery_task_logs.get(result.id, {}))


class AnnotationRunLeaseHeartbeatAbortTest(AnnotationRunLeaseTestCase):
    """ #1658: losing the lease must kill the in-flight subprocess, and holding it must not.

        The subprocess runs on a worker thread (as VEP does, blocking in communicate()) while _renew is
        driven from the test thread - all DB access stays on the connection inside TestCase's
        transaction, which a real heartbeat thread could not see. """

    def _run_subprocess_async(self, heartbeat, cmd):
        """ Start `cmd` under execute_cmd on a worker thread, registering it with the heartbeat exactly
            as dump_and_annotate_variants does. Returns (thread, result_holder) once the subprocess is
            registered and live. """
        registered = threading.Event()
        holder = {}

        def _register(process):
            heartbeat.set_process(process)
            if process is not None:
                registered.set()

        def _target():
            holder["result"] = execute_cmd(cmd, process_callback=_register)

        thread = threading.Thread(target=_target, daemon=True)
        thread.start()
        self.assertTrue(registered.wait(timeout=10), "subprocess was never handed to the heartbeat")
        return thread, holder

    def test_lost_lease_kills_running_subprocess(self):
        run = self._make_run(task_id="attempt-a")
        heartbeat = AnnotationRunLeaseHeartbeat(run, "attempt-a")
        thread, holder = self._run_subprocess_async(heartbeat, ["sleep", "60"])

        # Dispatcher reclaimed the run and handed it to another attempt.
        AnnotationRun.objects.filter(pk=run.pk).update(task_id="attempt-b")
        started = time.monotonic()
        heartbeat._renew()

        thread.join(timeout=10)
        self.assertFalse(thread.is_alive(), "subprocess outlived the lost lease")
        self.assertLess(time.monotonic() - started, 10, "abort did not happen promptly")
        self.assertNotEqual(holder["result"].return_code, 0,
                            "killed VEP must return non-zero so the pipeline fails this attempt")

    def test_held_lease_leaves_subprocess_alone(self):
        """ Counterpart: a healthy renew must never touch the subprocess. A false abort silently destroys
            a completed multi-hour VEP run and surfaces only as an unexplained pipeline failure. """
        run = self._make_run(task_id="attempt-a")
        original_expires = run.lease_expires
        heartbeat = AnnotationRunLeaseHeartbeat(run, "attempt-a")
        thread, holder = self._run_subprocess_async(heartbeat, ["sleep", "0.5"])

        heartbeat._renew()  # still ours - renews, must not kill

        thread.join(timeout=10)
        self.assertFalse(thread.is_alive())
        self.assertEqual(holder["result"].return_code, 0, "healthy run's subprocess was killed")

        run.refresh_from_db()
        self.assertGreater(run.lease_expires, original_expires, "lease was not renewed")


class AnnotationRunDumpPathUniquenessTest(AnnotationRunLeaseTestCase):
    """ #1658: unique per-task paths are the standalone correctness fix - each VEP owns its own file, so
        neither attempt can corrupt the other regardless of who wins the DB write race. """

    def test_distinct_task_tokens_give_distinct_dump_and_annotated_paths(self):
        run = self._make_run()
        dump_a = run.get_dump_filename(task_token="attempt-a")
        dump_b = run.get_dump_filename(task_token="attempt-b")
        self.assertNotEqual(dump_a, dump_b)

        # The annotated VCF name is derived from the dump name (dump_and_annotate_variants), so the
        # dump stems must differ, not merely the full paths.
        self.assertNotEqual(name_from_filename(dump_a), name_from_filename(dump_b))
        self.assertEqual(os.path.dirname(dump_a), os.path.dirname(dump_b))

    def test_pipeline_dumps_each_attempt_to_its_own_path(self):
        """ Guards the single call site that threads task_id into the dump path. Dropping that one kwarg
            silently restores the collision the issue is about, and nothing else would notice. """
        run = self._make_run(task_id="attempt-a")
        paths = []
        for task_id in ("attempt-a", "attempt-b"):
            run.task_id = task_id
            run.save()
            with mock.patch.object(annotate_variants_module, "_unannotated_variants_to_vcf", return_value=0), \
                    mock.patch.object(annotate_variants_module, "vep_check_command_line_version_match"):
                dump_and_annotate_variants(run)
            paths.append(run.vcf_dump_filename)

        self.assertNotEqual(paths[0], paths[1],
                            "two attempts on the same run dumped to the same path")

    def test_external_dump_path_is_stable_regardless_of_task_id(self):
        """ #1568 external dumps are copied to other machines and matched via their sidecar metadata, so
            their name must not acquire a local task token. """
        run = self._make_run(external=True)
        dump_before = run.get_dump_filename()

        run.task_id = "some-local-task"
        run.save()
        self.assertEqual(run.get_dump_filename(), dump_before,
                         "external dump path changed once the run held a task_id")

    def test_external_dump_and_metadata_sidecar_share_a_stem(self):
        """ get_dump_metadata_filename takes no task token while get_dump_filename does - they must still
            land next to each other for external matching to work. """
        run = self._make_run(external=True)
        dump = run.get_dump_filename()
        metadata = run.get_dump_metadata_filename()

        self.assertTrue(dump.endswith(".vcf.gz"))
        self.assertTrue(metadata.endswith(".meta.json"))
        self.assertEqual(dump[:-len(".vcf.gz")], metadata[:-len(".meta.json")],
                         "dump VCF and its metadata sidecar no longer share a stem")


class AnnotationRunTaskLogTest(AnnotationRunLeaseTestCase):
    """ celery_task_logs is keyed by task_id specifically so a reclaimed run keeps the logs of every
        attempt - which is the audit trail for diagnosing a reclaim in the first place. """

    def test_set_task_log_round_trips(self):
        run = self._make_run(task_id="attempt-a")
        run.set_task_log("start", "2026-07-21T00:00:00")
        run.save()

        run.refresh_from_db()
        self.assertEqual(run.celery_task_logs.get("attempt-a", {}).get("start"), "2026-07-21T00:00:00")

    def test_set_task_log_accepts_the_datetimes_the_pipeline_logs(self):
        """ Every real caller logs a timezone.now() (annotate_variants "start"/"end", import_annotation_run
            "import_start"/"import_end"). celery_task_logs is a plain JSONField with no encoder, so whatever
            set_task_log stores has to already be JSON-serializable. """
        run = self._make_run(task_id="attempt-a")
        run.set_task_log("start", timezone.now())
        run.save()

        run.refresh_from_db()
        self.assertIsNotNone(run.celery_task_logs.get("attempt-a", {}).get("start"))

    def test_set_task_log_keeps_entries_from_earlier_attempts(self):
        run = self._make_run(task_id="attempt-a")
        run.set_task_log("start", "first")
        run.save()

        run.task_id = "attempt-b"
        run.set_task_log("start", "second")
        run.save()

        run.refresh_from_db()
        self.assertEqual(run.celery_task_logs.get("attempt-a", {}).get("start"), "first",
                         "losing attempt's log was lost when the run was reclaimed")
        self.assertEqual(run.celery_task_logs.get("attempt-b", {}).get("start"), "second")
