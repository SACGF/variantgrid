import json
from contextlib import redirect_stdout
from importlib import import_module
from io import StringIO

from django.apps import apps
from django.core.management import CommandError, call_command
from django.test import TestCase

from manual import gates
from manual.models import ManualGateSatisfied, ManualMigrationAttempt, ManualMigrationOutstanding, \
    ManualMigrationRequired, ManualMigrationTask
from manual.operations.manual_operations import ManualOperation
from scripts.migrator.migrator import Migrator, run_scheduler

# Migration module names start with a digit, so they can't be imported with normal import syntax.
complete_obsolete_tasks = import_module(
    "manual.migrations.0004_complete_obsolete_manual_tasks").complete_obsolete_tasks
backfill_requires = import_module(
    "manual.migrations.0005_backfill_manual_task_requires").backfill_requires


class GateResolutionTest(TestCase):
    """ How a task's requires resolve to blocked/satisfied - the logic gating auto-run. """

    def _request(self, task):
        ManualMigrationRequired.objects.create(task=task)

    def test_manual_gate_blocks_until_operator_confirms(self):
        gate = "variant-annotation-current"  # a ManualGate - only satisfied by operator confirmation
        self.assertEqual(gates.blocked_by([gate]), [gate])
        ManualGateSatisfied.objects.create(name=gate)
        self.assertEqual(gates.blocked_by([gate]), [])

    def test_unknown_gate_blocks_fail_safe(self):
        # A typo'd/unregistered gate must block, not silently pass (never run something ungated).
        self.assertEqual(gates.blocked_by(["typo-gate"]), ["typo-gate"])

    def test_after_task_gate_tracks_prerequisite(self):
        dep = ManualMigrationTask.objects.create(id="manage*prereq_step")
        after = f"{gates.AFTER_PREFIX}{dep.id}"

        self.assertEqual(gates.blocked_by([after]), [])          # never requested -> nothing to wait on
        self._request(dep)
        self.assertEqual(gates.blocked_by([after]), [after])     # outstanding -> blocks
        ManualMigrationAttempt.objects.create(task=dep, requires_retry=False)
        self.assertEqual(gates.blocked_by([after]), [])          # completed -> released

    def test_operation_persists_requires_to_task(self):
        ManualOperation.operation_manage(["match_patient_phenotypes", "--clear"],
                                         requires=["ontology-imported"]).run(apps)
        task = ManualMigrationTask.objects.get(pk="manage*match_patient_phenotypes --clear")
        self.assertEqual(task.requires, ["ontology-imported"])

    def test_manual_gate_command_only_confirms_manual_gates(self):
        with self.assertRaises(CommandError):
            call_command("manual_gate", "--satisfy", "ontology-imported")  # auto gate can't be confirmed
        with self.assertRaises(CommandError):
            call_command("manual_gate", "--satisfy", "not-a-gate")
        call_command("manual_gate", "--satisfy", "variant-annotation-current")
        self.assertTrue(ManualGateSatisfied.is_satisfied("variant-annotation-current"))


class OutstandingRunnableTest(TestCase):
    """ manual_outstanding's per-task runnable/blocked_by/command_exists - what the migrator acts on. """

    def _request(self, task):
        ManualMigrationRequired.objects.create(task=task)

    def _outstanding_by_id(self):
        out = StringIO()
        with redirect_stdout(out):  # command print()s its JSON (migrator reads it via Popen)
            call_command("manual_outstanding")
        return {t["id"]: t for t in json.loads(out.getvalue())["tasks"]}

    def test_runnable_only_when_manage_command_exists_and_ungated(self):
        ready = ManualMigrationTask.objects.create(id="manage*migrate")  # real command, no gate
        gated = ManualMigrationTask.objects.create(
            id="manage*showmigrations", requires=["variant-annotation-current"])
        missing = ManualMigrationTask.objects.create(id="manage*deleted_command")  # command gone
        human = ManualMigrationTask.objects.create(id="other*do a manual thing")
        for t in (ready, gated, missing, human):
            self._request(t)

        by_id = self._outstanding_by_id()

        self.assertTrue(by_id[ready.id]["runnable"])
        self.assertFalse(by_id[gated.id]["runnable"])                                   # blocked by gate
        self.assertEqual(by_id[gated.id]["blocked_by"], ["variant-annotation-current"])
        self.assertFalse(by_id[missing.id]["command_exists"])
        self.assertFalse(by_id[missing.id]["runnable"])                                 # command missing
        self.assertFalse(by_id[human.id]["runnable"])                                   # 'other' never runs
        self.assertIsNone(by_id[human.id]["command_exists"])                            # no command -> not obsolete

    def test_menu_status_line_flags_only_missing_manage_commands(self):
        # Regression: 'other' human steps have command_exists=None and must NOT read as "obsolete command".
        missing = Migrator.subcommand_for_json(
            {"id": "manage*deleted_cmd", "category": "manage", "line": "deleted_cmd",
             "command_exists": False, "blocked_by": []})
        human = Migrator.subcommand_for_json(
            {"id": 'other*"do a thing"', "category": "other", "line": '"do a thing"',
             "command_exists": None, "blocked_by": []})
        self.assertIn("obsolete", missing.status_line())
        self.assertIsNone(human.status_line())


class ObsoleteCleanupTest(TestCase):
    """ complete_obsolete_tasks (manual/0004): retire dead tasks, never touch ones we keep. """

    def _request(self, task):
        ManualMigrationRequired.objects.create(task=task)

    def _outstanding(self, task):
        return ManualMigrationOutstanding.outstanding_task(task) is not None

    def test_completes_obsolete_but_leaves_kept_tasks(self):
        # command deleted -> completed generically (can't run anyway)
        missing_cmd = ManualMigrationTask.objects.create(id="manage*one_off_calc_variant_end")
        # command still exists, but this exact variant is a vetted skip
        skip_variant = ManualMigrationTask.objects.create(id="manage*gene_annotation --add-missing-omim")
        # free-text reminder retired by exact id
        reminder = ManualMigrationTask.objects.create(
            id='other*"Import dbNSFP gene annotation (see annotation page)"')
        # kept: a sibling of the skipped variant, a live command, and an unrelated human step
        kept_sibling = ManualMigrationTask.objects.create(id="manage*gene_annotation --new-releases")
        kept_command = ManualMigrationTask.objects.create(id="manage*migrate")
        kept_reminder = ManualMigrationTask.objects.create(id='other*"A human step we keep"')
        for t in (missing_cmd, skip_variant, reminder, kept_sibling, kept_command, kept_reminder):
            self._request(t)

        complete_obsolete_tasks(apps, None)

        self.assertFalse(self._outstanding(missing_cmd))
        self.assertFalse(self._outstanding(skip_variant))
        self.assertFalse(self._outstanding(reminder))
        self.assertTrue(self._outstanding(kept_sibling))    # exact-id match must not hit siblings
        self.assertTrue(self._outstanding(kept_command))
        self.assertTrue(self._outstanding(kept_reminder))


class SchedulerTest(TestCase):
    """ The auto-manage scheduler: run unblocked tasks, re-evaluate between passes, stop on failure. """

    def _request(self, task):
        ManualMigrationRequired.objects.create(task=task)

    def test_run_scheduler_reevaluates_between_passes(self):
        # B only becomes runnable once A is done - proves it re-fetches rather than snapshotting once.
        done, graph = set(), {"A": set(), "B": {"A"}, "C": {"B"}}
        ran = []

        def fetch():
            return [t for t in ("A", "B", "C") if t not in done and graph[t] <= done]

        def run_ok(task):
            ran.append(task)
            done.add(task)
            return True

        completed, failed = run_scheduler(fetch, run_ok)
        self.assertEqual(ran, ["A", "B", "C"])
        self.assertIsNone(failed)

    def test_run_scheduler_stops_on_first_failure(self):
        attempted = []

        def fetch():
            return ["X", "Y"] if "X" not in attempted else []

        def run(task):
            attempted.append(task)
            return task != "X"  # X fails

        completed, failed = run_scheduler(fetch, run)
        self.assertEqual(attempted, ["X"])  # Y in the same pass is not attempted after X fails
        self.assertEqual(failed, "X")

    def test_after_chain_runs_in_insertion_order_end_to_end(self):
        # Backfilled after-gates (see 0005) must make the fix_variant_matching-style chain run strictly
        # in order. Uses the real run_scheduler + manual_outstanding + gate resolution; no gates besides
        # the after-chain, so nothing external blocks it.
        a = ManualMigrationTask.objects.create(id="manage*migrate")
        b = ManualMigrationTask.objects.create(
            id="manage*showmigrations", requires=["after:manage*migrate"])
        c = ManualMigrationTask.objects.create(
            id="manage*makemigrations", requires=["after:manage*showmigrations"])
        for t in (a, b, c):
            self._request(t)
        chain_ids = {a.id, b.id, c.id}  # the migrated test DB has other outstanding tasks - ignore them

        def fetch_runnable():
            out = StringIO()
            with redirect_stdout(out):
                call_command("manual_outstanding")
            return [t for t in json.loads(out.getvalue())["tasks"]
                    if t["runnable"] and t["id"] in chain_ids]

        ran = []

        def run_task(task):
            ran.append(task["line"])
            ManualMigrationAttempt.objects.create(
                task=ManualMigrationTask.objects.get(pk=task["id"]), requires_retry=False)
            return True

        run_scheduler(fetch_runnable, run_task)
        self.assertEqual(ran, ["migrate", "showmigrations", "makemigrations"])
