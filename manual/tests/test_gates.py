import json
from contextlib import redirect_stdout
from io import StringIO

from django.apps import apps
from django.core.management import CommandError, call_command
from django.test import TestCase

from importlib import import_module

from manual import gates
from manual.models import ManualGateSatisfied, ManualMigrationAttempt, ManualMigrationOutstanding, \
    ManualMigrationRequired, ManualMigrationTask
from manual.operations.manual_operations import ManualOperation

# Migration module names start with a digit, so they can't be imported with normal import syntax.
complete_obsolete_tasks = import_module(
    "manual.migrations.0004_complete_obsolete_manual_tasks").complete_obsolete_tasks
backfill_requires = import_module(
    "manual.migrations.0005_backfill_manual_task_requires").backfill_requires


class ManualGateTest(TestCase):

    def _request(self, task: ManualMigrationTask):
        """ Mark a task as outstanding (a migration would create this during migrate). """
        ManualMigrationRequired.objects.create(task=task)

    def test_manual_gate_blocks_until_satisfied(self):
        gate = "variant-annotation-current"  # a ManualGate
        self.assertEqual(gates.blocked_by([gate]), [gate])
        ManualGateSatisfied.objects.create(name=gate)
        self.assertEqual(gates.blocked_by([gate]), [])

    def test_unknown_gate_always_blocks(self):
        self.assertEqual(gates.blocked_by(["typo-gate"]), ["typo-gate"])

    def test_after_task_gate(self):
        dep = ManualMigrationTask.objects.create(id="manage*test_gate_dep")
        after = f"{gates.AFTER_PREFIX}{dep.id}"

        # Never requested -> nothing to wait for -> satisfied
        self.assertEqual(gates.blocked_by([after]), [])

        # Requested but not done -> blocks
        self._request(dep)
        self.assertEqual(gates.blocked_by([after]), [after])

        # A recorded success clears it (no longer outstanding)
        ManualMigrationAttempt.objects.create(task=dep, requires_retry=False)
        self.assertEqual(gates.blocked_by([after]), [])

    def test_operation_persists_requires(self):
        op = ManualOperation.operation_manage(["match_patient_phenotypes", "--clear"],
                                              requires=["ontology-imported"])
        op.run(apps)  # run(apps).get_model works against the live registry in tests
        task = ManualMigrationTask.objects.get(pk="manage*match_patient_phenotypes --clear")
        self.assertEqual(task.requires, ["ontology-imported"])

    def test_outstanding_reports_blocked_and_runnable(self):
        blocked_task = ManualMigrationTask.objects.create(
            id="manage*test_gate_calc_stats", requires=["variant-annotation-current"])
        self._request(blocked_task)
        ready_task = ManualMigrationTask.objects.create(id="manage*test_gate_ready")
        self._request(ready_task)

        outstanding = {o.task.id: o.to_json() for o in ManualMigrationOutstanding.outstanding_tasks()}
        blocked_json = outstanding[blocked_task.id]
        blocked_json["blocked_by"] = gates.blocked_by(blocked_json["requires"])
        self.assertEqual(blocked_json["blocked_by"], ["variant-annotation-current"])

        ready_json = outstanding[ready_task.id]
        self.assertEqual(gates.blocked_by(ready_json["requires"]), [])

    def test_manual_gate_command_satisfy(self):
        with self.assertRaises(CommandError):
            call_command("manual_gate", "--satisfy", "ontology-imported")  # auto gate can't be confirmed
        with self.assertRaises(CommandError):
            call_command("manual_gate", "--satisfy", "not-a-gate")

        call_command("manual_gate", "--satisfy", "variant-annotation-current", "--note", "done")
        self.assertTrue(ManualGateSatisfied.is_satisfied("variant-annotation-current"))

    def test_manual_outstanding_command_marks_runnable(self):
        blocked = ManualMigrationTask.objects.create(
            id="manage*migrate", requires=["variant-annotation-current"])  # 'migrate' is a real command
        self._request(blocked)
        ready = ManualMigrationTask.objects.create(id="manage*showmigrations")  # real command, no gate
        self._request(ready)
        other = ManualMigrationTask.objects.create(id="other*test_gate_manual_step")
        self._request(other)

        out = StringIO()
        with redirect_stdout(out):  # manual_outstanding print()s its JSON (migrator reads it via Popen)
            call_command("manual_outstanding")
        by_id = {t["id"]: t for t in json.loads(out.getvalue())["tasks"]}

        self.assertFalse(by_id[blocked.id]["runnable"])
        self.assertEqual(by_id[blocked.id]["blocked_by"], ["variant-annotation-current"])
        self.assertTrue(by_id[ready.id]["runnable"])
        self.assertFalse(by_id[other.id]["runnable"])  # 'other' category never auto-runs

    def test_missing_command_is_not_runnable(self):
        # 'fix_something' isn't a real management command -> obsolete/orphaned row, never auto-run
        orphan = ManualMigrationTask.objects.create(id="manage*test_gate_no_such_command")
        self._request(orphan)
        real = ManualMigrationTask.objects.create(id="manage*migrate")  # 'migrate' is a real command
        self._request(real)

        out = StringIO()
        with redirect_stdout(out):
            call_command("manual_outstanding")
        by_id = {t["id"]: t for t in json.loads(out.getvalue())["tasks"]}

        self.assertFalse(by_id[orphan.id]["command_exists"])
        self.assertFalse(by_id[orphan.id]["runnable"])
        self.assertTrue(by_id[real.id]["command_exists"])
        self.assertTrue(by_id[real.id]["runnable"])

    def test_cleanup_migration_completes_obsolete_tasks(self):
        obsolete = ManualMigrationTask.objects.create(id="manage*one_off_calc_variant_end")
        self._request(obsolete)
        keep = ManualMigrationTask.objects.create(id="manage*migrate")  # not obsolete
        self._request(keep)
        self.assertIsNotNone(ManualMigrationOutstanding.outstanding_task(obsolete))

        complete_obsolete_tasks(apps, None)

        self.assertIsNone(ManualMigrationOutstanding.outstanding_task(obsolete))  # cleared
        self.assertIsNotNone(ManualMigrationOutstanding.outstanding_task(keep))   # untouched

    def test_cleanup_migration_completes_obsolete_task_id_but_not_kept_steps(self):
        # command still exists -> matched by exact id, must NOT clear commands/variants we keep
        omim = ManualMigrationTask.objects.create(id="manage*gene_annotation --add-missing-omim")
        self._request(omim)
        kept_variant = ManualMigrationTask.objects.create(id="manage*gene_annotation --new-releases")
        self._request(kept_variant)
        link_transcripts = ManualMigrationTask.objects.create(id="manage*fix_annotation_link_transcripts")
        self._request(link_transcripts)

        complete_obsolete_tasks(apps, None)

        self.assertIsNone(ManualMigrationOutstanding.outstanding_task(omim))            # cleared
        self.assertIsNotNone(ManualMigrationOutstanding.outstanding_task(kept_variant))     # untouched
        self.assertIsNotNone(ManualMigrationOutstanding.outstanding_task(link_transcripts))  # kept

    def test_backfill_requires_migration(self):
        extra = ManualMigrationTask.objects.create(id="manage*fix_variant_matching --extra")
        reval = ManualMigrationTask.objects.create(id="manage*fix_variant_matching --revalidate_chgvs")
        phenos = ManualMigrationTask.objects.create(id="manage*match_patient_phenotypes --clear")

        backfill_requires(apps, None)

        for t in (extra, reval, phenos):
            t.refresh_from_db()
        self.assertEqual(extra.requires, ["cdot-current"])
        # Chained in insertion order: revalidate_chgvs runs after --extra
        self.assertEqual(reval.requires,
                         ["cdot-current", "after:manage*fix_variant_matching --extra"])
        self.assertEqual(phenos.requires, ["ontology-imported"])

        # The after-gate blocks revalidate_chgvs while --extra is still outstanding, and clears once done
        self._request(extra)
        self.assertIn("after:manage*fix_variant_matching --extra", gates.blocked_by(reval.requires))
        ManualMigrationAttempt.objects.create(task=extra, requires_retry=False)
        self.assertNotIn("after:manage*fix_variant_matching --extra", gates.blocked_by(reval.requires))
