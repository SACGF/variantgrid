from contextlib import redirect_stdout
from importlib import import_module
from io import StringIO
from unittest.mock import patch

from django.apps import apps
from django.core.management import CommandError, call_command
from django.test import TestCase

from manual import gates
from manual.models import (
    ManualGateSatisfied,
    ManualMigrationAttempt,
    ManualMigrationOutstanding,
    ManualMigrationRequired,
    ManualMigrationTask,
)
from manual.operations.manual_operations import ManualOperation
from manual.upgrader import (
    ObsoleteStep,
    StepStatus,
    Upgrader,
    outstanding_tasks_json,
    record_attempt,
    run_scheduler,
)

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
    """ outstanding_tasks_json's per-task runnable/blocked_by/command_exists - what the upgrader acts on. """

    def _request(self, task):
        ManualMigrationRequired.objects.create(task=task)

    def _outstanding_by_id(self):
        return {t["id"]: t for t in outstanding_tasks_json()}

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
        missing = Upgrader.step_for_task(
            {"id": "manage*deleted_cmd", "category": "manage", "line": "deleted_cmd",
             "command_exists": False, "blocked_by": []})
        human = Upgrader.step_for_task(
            {"id": 'other*"do a thing"', "category": "other", "line": '"do a thing"',
             "command_exists": None, "blocked_by": []})
        self.assertEqual(missing.status_tag(), "[OBSOLETE]")
        self.assertIsNone(human.status_tag())
        self.assertIsNone(human.status_line())

    def test_menu_tags_blocked_task_on_its_own_line(self):
        # The tag has to sit on the task's own menu line - an indented detail line alone reads as
        # belonging to whichever task is printed next.
        blocked = Upgrader.step_for_task(
            {"id": "manage*calculate_sample_stats", "category": "manage", "line": "calculate_sample_stats",
             "command_exists": True, "blocked_by": ["variant-annotation-current"]})
        self.assertEqual(blocked.status_tag(), "[BLOCKED]")
        self.assertIn("variant-annotation-current", blocked.status_line())

    def test_selecting_obsolete_task_offers_to_mark_it_complete(self):
        # Selecting an obsolete step used to shell out to a command that no longer exists and fail,
        # with no way to retire it from the upgrader.
        task = ManualMigrationTask.objects.create(id="manage*deleted_cmd")
        self._request(task)
        obsolete = Upgrader.step_for_task(
            {"id": task.id, "category": "manage", "line": "deleted_cmd",
             "command_exists": False, "blocked_by": []})
        self.assertIsInstance(obsolete, ObsoleteStep)

        out = StringIO()
        with redirect_stdout(out), patch("builtins.input", return_value="y"):
            result = obsolete.run()
        self.assertEqual(result.status, StepStatus.SUCCESS)
        self.assertTrue(result.note)

        # what the upgrader does with that success - records it, which retires the task
        record_attempt(task.id, note=result.note)
        self.assertIsNone(ManualMigrationOutstanding.outstanding_task(task))

    def test_backing_out_of_an_obsolete_task_leaves_it_outstanding(self):
        task = ManualMigrationTask.objects.create(id="manage*deleted_cmd")
        self._request(task)
        obsolete = ObsoleteStep("deleted_cmd").using(task_id=task.id)

        out = StringIO()
        with redirect_stdout(out), patch("builtins.input", return_value="x"):
            result = obsolete.run()
        self.assertEqual(result.status, StepStatus.SKIP)  # SKIP -> upgrader records no attempt
        self.assertIsNotNone(ManualMigrationOutstanding.outstanding_task(task))


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
    """ The auto-manage scheduler: run unblocked tasks, re-evaluate between passes, carry on past failure. """

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

        _, failed = run_scheduler(fetch, run_ok)
        self.assertEqual(ran, ["A", "B", "C"])
        self.assertEqual(failed, [])

    def test_run_scheduler_carries_on_past_failure_without_retrying(self):
        # X fails and stays outstanding (so fetch keeps offering it); Z is gated after X so never offered
        done, graph = set(), {"X": set(), "Y": set(), "Z": {"X"}}
        attempted = []

        def fetch():
            return [t for t in ("X", "Y", "Z") if t not in done and graph[t] <= done]

        def run(task):
            attempted.append(task)
            if task == "X":
                return False
            done.add(task)
            return True

        _, failed = run_scheduler(fetch, run)
        self.assertEqual(attempted, ["X", "Y"])
        self.assertEqual(failed, ["X"])

    def test_after_chain_runs_in_insertion_order_end_to_end(self):
        # Backfilled after-gates (see 0005) must make the fix_variant_matching-style chain run strictly
        # in order. Uses the real run_scheduler + outstanding_tasks_json + gate resolution; no gates besides
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
            return [t for t in outstanding_tasks_json() if t["runnable"] and t["id"] in chain_ids]

        ran = []

        def run_task(task):
            ran.append(task["line"])
            ManualMigrationAttempt.objects.create(
                task=ManualMigrationTask.objects.get(pk=task["id"]), requires_retry=False)
            return True

        run_scheduler(fetch_runnable, run_task, key=lambda task: task["id"])
        self.assertEqual(ran, ["migrate", "showmigrations", "makemigrations"])
