import sys
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

from django.test import TestCase

from manual.models import ManualMigrationAttempt, ManualMigrationRequired, ManualMigrationTask
from manual.upgrader import (
    PythonStep,
    StepResult,
    StepStatus,
    Upgrader,
    parse_selection,
)


def _step(task_id: str, method, requires=None) -> PythonStep:
    step = PythonStep(method, task_id).using(task_id=task_id)
    step.requires = requires
    return step


def _fail():
    raise ValueError("boom")


class ParseSelectionTest(TestCase):

    def test_keys_and_ranges_in_order_given(self):
        migrations = [PythonStep(StepResult.success, key).using(key=key)
                      for key in ("m", "c", "1", "2", "3", "4")]
        selected = parse_selection("c, 2-3 m", migrations)
        self.assertEqual([m.key for m in selected], ["c", "2", "3", "m"])
        with self.assertRaises(ValueError):
            parse_selection("2-5", migrations)  # 5 isn't on the menu


class RunSelectionTest(TestCase):

    def _run(self, upgrader: Upgrader, migrations, **kwargs):
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            return upgrader.run_selection(migrations, then="exit", **kwargs)

    def _attempts(self, task_id: str) -> list[bool]:
        return list(ManualMigrationAttempt.objects.filter(task_id=task_id).values_list("requires_retry", flat=True))

    def test_carries_on_past_failure_skipping_what_is_gated_after_it(self):
        for task_id in ("manage*a", "manage*b", "manage*c"):
            ManualMigrationRequired.objects.create(task=ManualMigrationTask.objects.create(id=task_id))
        ran = []
        a = _step("manage*a", _fail)
        b = _step("manage*b", lambda: ran.append("b") or StepResult.success(), requires=["after:manage*a"])
        c = _step("manage*c", lambda: ran.append("c") or StepResult.success())

        failed = self._run(Upgrader(), [a, b, c])

        self.assertEqual(failed, [a])
        self.assertEqual(ran, ["c"])
        self.assertEqual(self._attempts("manage*a"), [True])  # recorded as needing a retry
        self.assertEqual(self._attempts("manage*b"), [])      # skipped, not attempted
        self.assertEqual(self._attempts("manage*c"), [False])

    def test_stop_on_failure(self):
        ran = []
        failed = self._run(Upgrader(), [_step("manage*a", _fail),
                                        _step("manage*c", lambda: ran.append("c") or StepResult.success())],
                           stop_on_failure=True)
        self.assertEqual(len(failed), 1)
        self.assertEqual(ran, [])

    def test_command_exiting_zero_is_a_success(self):
        upgrader = Upgrader()
        with redirect_stdout(StringIO()):
            ok = upgrader.run_step(_step("manage*a", lambda: sys.exit(0)))
            bad = upgrader.run_step(_step("manage*b", lambda: sys.exit(2)))
        self.assertEqual(ok.status, StepStatus.SUCCESS)
        self.assertEqual(bad.status, StepStatus.FAILURE)
