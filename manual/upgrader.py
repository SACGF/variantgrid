"""
The deploy upgrader behind `scripts/upgrade.sh` (`manage.py upgrade`): one Django process runs the standard
steps (git pull, migrate, collectstatic, deployment_check, deployed) and the outstanding ManualMigrationTasks
in-process with call_command, records each attempt, and offers a menu for the rest. A pull that moves HEAD
re-execs the process (`Upgrader.restart_after_pull`) so the remaining steps run on the new code.

Entry points: Upgrader (quick / auto_manage / run_selection / prompt), outstanding_tasks_json (also behind
`manage.py manual_outstanding`), record_attempt and run_scheduler.
"""
import json
import os
import re
import shlex
import socket
import subprocess
import sys
import traceback
from collections.abc import Callable
from enum import Enum, auto
from functools import partial
from typing import Any, Optional

import requests
from django.conf import settings
from django.core.management import call_command, get_commands
from django.db import connections

from library.constants import MINUTE_SECS
from library.git import Git
from manual.gates import AFTER_PREFIX, blocked_by
from manual.models import ManualMigrationAttempt, ManualMigrationOutstanding, ManualMigrationTask
from variantgrid.settings.components.secret_settings import get_secret

YELLOW = "\033[93m"
COLOR_END = "\033[00m"


def print_color(color: str, skk: str):
    print((color + "{}" + COLOR_END).format(skk))


def color(the_color: str, text: str) -> str:
    return the_color + text + COLOR_END


print_red = partial(print_color, "\033[91m")
print_green = partial(print_color, "\033[92m")
print_yellow = partial(print_color, YELLOW)
print_light_purple = partial(print_color, "\033[94m")
print_purple = partial(print_color, "\033[95m")
print_cyan = partial(print_color, "\033[96m")


def print_restart_reminder():
    print_yellow("If this upgrade pulled new code, restart the services so they pick it up:")
    print_yellow("    ctrl-d (back to an admin user)")
    print_yellow("    sudo ./scripts/restart_services.sh")


def outstanding_tasks_json() -> list[dict[str, Any]]:
    """ Each outstanding task's to_json() plus 'command_exists', 'blocked_by' and 'runnable' """
    known_commands = set(get_commands())
    task_list = []
    for outstanding_task in ManualMigrationOutstanding.outstanding_tasks():
        task_json = outstanding_task.to_json()
        is_manage = task_json["category"] == "manage"
        # command_exists only applies to 'manage' tasks (None for 'other' human steps, which have no
        # command). A manage task whose command was deleted is an obsolete/orphaned row - never
        # auto-run it (it would crash), it's reported for review instead.
        command_exists = task_json["line"].split()[0] in known_commands if is_manage else None
        task_json["command_exists"] = command_exists
        blocking = blocked_by(task_json.get("requires"))
        task_json["blocked_by"] = blocking
        # Only 'manage' tasks are auto-runnable, and only once nothing blocks them and the command exists.
        task_json["runnable"] = bool(is_manage and command_exists and not blocking)
        task_list.append(task_json)
    return task_list


def record_attempt(task_id: str, success: bool = True, note: Optional[str] = None,
                   version: Optional[str] = None):
    task, _ = ManualMigrationTask.objects.get_or_create(id=task_id)
    ManualMigrationAttempt.objects.create(task=task, note=note, source_version=version,
                                          requires_retry=not success)


class MigrationStatus(Enum):
    SUCCESS = auto()
    FAILURE = auto()
    SKIP = auto()


class MigrationResult:

    def __init__(self, status: MigrationStatus, note: Optional[str] = None, new_code: bool = False):
        self.status = status
        self.note = note
        self.new_code = new_code  # a git pull moved HEAD - the rest has to run in a fresh process

    @staticmethod
    def success(note: Optional[str] = None):
        return MigrationResult(status=MigrationStatus.SUCCESS, note=note)

    @staticmethod
    def failure(note: Optional[str] = None):
        return MigrationResult(status=MigrationStatus.FAILURE, note=note)

    @staticmethod
    def skip():
        return MigrationResult(status=MigrationStatus.SKIP)


class SubMigration:

    def __init__(self):
        self.key = None
        self.task_id = None
        self.notes = None
        self.requires = None         # gate names (see manual.gates)
        self.blocked_by = None       # gate names not yet satisfied
        self.command_exists = None   # False for an obsolete 'manage' task whose command was deleted

    def using(self, key: Optional[str] = None, task_id: Optional[str] = None, notes: Optional[list[str]] = None):
        if key:
            self.key = key
        if task_id:
            self.task_id = task_id
        if notes:
            self.notes = notes
        return self

    @property
    def ref(self) -> str:
        """ Stable across menu refreshes (and a restart), unlike the numbered key of a manual task """
        return self.task_id or self.key

    def status_tag(self) -> Optional[str]:
        """ Short prefix for the menu entry itself, so the status can't be mistaken for the
            neighbouring task's. """
        if self.command_exists is False:
            return "[OBSOLETE]"
        if self.blocked_by:
            return "[BLOCKED]"
        return None

    def status_line(self) -> Optional[str]:
        """ Detail for the status_tag, printed indented under the menu entry, or None. """
        if self.command_exists is False:
            return "command no longer exists (selecting it here marks it complete)"
        if self.blocked_by:
            return f"waiting on: {', '.join(self.blocked_by)} (selecting it here still runs it)"
        return None

    def run(self) -> MigrationResult:
        return MigrationResult.skip()


class PythonSubMigration(SubMigration):

    def __init__(self, the_method: Callable[[], MigrationResult], description: str):
        super().__init__()
        self.the_method = the_method
        self.description = description

    def run(self):
        return self.the_method()

    def __str__(self):
        return self.description


class ManualSubMigration(SubMigration):

    def __init__(self, text: str):
        super().__init__()
        self.text = text

    def __str__(self):
        return "Manually " + self.text

    def run(self):
        while True:
            print(self)
            print("y: record success")
            print("n: record failure")
            print("x: back")
            selection = input("\033[95mPlease enter a selection: \033[00m")
            selection = selection.strip().lower()
            if selection == "y":
                return MigrationResult.success()
            if selection == "n":
                return MigrationResult.failure()
            if selection == "x":
                return MigrationResult.skip()
            print(f"Unexpected input - \"{selection}\"")


class ObsoleteSubMigration(SubMigration):
    """ A 'manage' task whose command no longer exists - running it can only fail, so selecting it
        offers to record it as complete (which takes it off the outstanding list) instead. """

    def __init__(self, line: str):
        super().__init__()
        self.line = line
        self.command_exists = False

    def __str__(self):
        return "python3 manage.py " + self.line

    def run(self) -> MigrationResult:
        while True:
            print_yellow(f"'{self.line}' no longer exists in this codebase, so it can't be run.")
            print("y: mark as complete (no longer required)")
            print("x: back")
            selection = input("\033[95mPlease enter a selection: \033[00m")
            selection = selection.strip().lower()
            if selection == "y":
                return MigrationResult.success(note="Obsolete - command no longer exists, marked complete in upgrader")
            if selection == "x":
                return MigrationResult.skip()
            print(f"Unexpected input - \"{selection}\"")


class GitPullSubMigration(SubMigration):
    """ Pulls, and if that moved HEAD installs the requirements it brought (upgrade.sh already installed
        them for the code we started on). The new code is loaded by restarting the process. """

    def __str__(self):
        return "git pull + install requirements"

    def run(self) -> MigrationResult:
        print_cyan(str(self))
        before = Git(settings.BASE_DIR).hash
        if (process := subprocess.run(["git", "pull"], cwd=settings.BASE_DIR, check=False)).returncode != 0:
            return MigrationResult.failure(f"'git pull' failed with error code {process.returncode}")
        if Git(settings.BASE_DIR).hash == before:
            return MigrationResult.success()

        install = os.path.join(settings.BASE_DIR, "scripts", "install_requirements.sh")
        if (process := subprocess.run([install], cwd=settings.BASE_DIR, check=False)).returncode != 0:
            return MigrationResult.failure(f"'{install}' failed with error code {process.returncode}")
        return MigrationResult(status=MigrationStatus.SUCCESS, new_code=True)


class ManageSubMigration(SubMigration):
    """ A management command, run in this process """

    def __init__(self, args: list[str]):
        super().__init__()
        self.args = args

    def __str__(self):
        return "python3 manage.py " + shlex.join(self.args)

    def run(self) -> MigrationResult:
        print_cyan(str(self))
        print_purple("-----------")
        try:
            call_command(*self.args)
        finally:
            print_purple("-----------")
        return MigrationResult.success()


def run_scheduler(fetch_runnable, run_task, key: Callable = lambda task: task, max_passes: int = 1000):
    """ Pure scheduler loop, decoupled from Django for testability.

        Repeatedly asks fetch_runnable() for the currently-runnable tasks and runs each, re-fetching
        after every pass so ordering gates (after:<task_id>) release as their prerequisites complete.
        Carries on past a failure: a failed task stays outstanding, so whatever is gated after it stays
        blocked, and it isn't retried this run. max_passes is a runaway guard.

        fetch_runnable() -> list[task];  run_task(task) -> bool (True on success).
        Returns (ran, failed): tasks run in order, and those that failed. """
    ran = []
    failed = []
    attempted = set()
    for _ in range(max_passes):
        runnable = [task for task in fetch_runnable() if key(task) not in attempted]
        if not runnable:
            break
        for task in runnable:
            attempted.add(key(task))
            ran.append(task)
            if not run_task(task):
                failed.append(task)
    return ran, failed


def parse_selection(selection: str, migrations: list[SubMigration]) -> list[SubMigration]:
    """ 'm, c 2-5' -> those migrations, in the order given. Raises ValueError naming anything unknown """
    by_key = {migration.key: migration for migration in migrations}
    selected = []
    for token in re.split(r"[\s,]+", selection.strip()):
        if not token:
            continue
        if m := re.fullmatch(r"(\d+)-(\d+)", token):
            keys = [str(i) for i in range(int(m.group(1)), int(m.group(2)) + 1)]
        else:
            keys = [token]
        for key in keys:
            if key not in by_key:
                raise ValueError(f"Unknown step '{key}'")
            selected.append(by_key[key])
    return selected


class Upgrader:
    STANDARD_KEYS = "g, m, r, c, k, d"

    def __init__(self):
        self.git_version = Git(settings.BASE_DIR).hash

    def standard_migrations(self) -> list[SubMigration]:
        return [
            GitPullSubMigration().using(key="g", task_id="git*pull"),
            ManageSubMigration(["migrate"]).using(key="m", task_id="manage*migrate"),
            ManageSubMigration(["collectstatic_js_reverse"]).using(key="r", task_id="manage*collectstatic_js_reverse"),
            # collectstatic without warning for conflicting files has been an issue for 6 years
            # see https://code.djangoproject.com/ticket/26583 maybe it'll get fixed soon? For now (since we've never had
            # a problem) just turn off all verbosity
            # --clear so a moved file whose mtime looks unchanged doesn't leave a stale copy behind. The wrapper also
            # drops the compressor's cached tags, which point at the bundles --clear just deleted
            ManageSubMigration(["collectstatic_clean_compressor", "-v", "0", "--noinput", "--clear"]).using(
                key="c", task_id="manage*collectstatic_clean_compressor"),
            ManageSubMigration(["deployment_check", "--die-if-invalid", "--quiet"]).using(
                key="k", task_id="manage*deployment_check"),
            PythonSubMigration(self.notify_deployed, "record deployment (Rollbar + manage.py deployed)").using(key="d"),
        ]

    @staticmethod
    def subcommand_for_json(task: dict) -> SubMigration:
        task_id = task["id"]
        line = task["line"]
        notes = task.get("notes")
        if task["category"] == "manage":
            if task.get("command_exists") is False:
                sub_migration = ObsoleteSubMigration(line)
            else:
                sub_migration = ManageSubMigration(shlex.split(line))
        else:
            sub_migration = ManualSubMigration(line)
        sub_migration.using(task_id=task_id, notes=notes)
        sub_migration.requires = task.get("requires")
        sub_migration.blocked_by = task.get("blocked_by")
        sub_migration.command_exists = task.get("command_exists")
        return sub_migration

    def custom_migrations(self) -> list[SubMigration]:
        return [Upgrader.subcommand_for_json(task).using(key=str(i))
                for i, task in enumerate(outstanding_tasks_json(), start=1)]

    def menu_migrations(self) -> list[SubMigration]:
        return self.standard_migrations() + self.custom_migrations()

    def notify_deployed(self) -> MigrationResult:
        rollbar_token = get_secret("ROLLBAR.access_token", mandatory=False)
        if not rollbar_token:
            print_red("No rollbar token found")
            return MigrationResult.failure()

        data = {
            "access_token": rollbar_token,
            "environment": socket.gethostname().lower().split('.')[0].replace('-', ''),
            "revision": self.git_version,
            "local_username": os.getenv("USER", ""),
        }
        try:
            response = requests.post("https://api.rollbar.com/api/1/deploy/", data=data, timeout=MINUTE_SECS)
            if response.status_code == 200:
                call_command("deployed")
            else:
                print(f"Failed to record deployment in Rollbar. Response: {response.text}")
        except requests.RequestException as e:
            print(f"Error recording deployment in Rollbar: {e!s}")
        return MigrationResult.success()

    def run_step(self, migration: SubMigration) -> MigrationResult:
        """ Runs one step, records the attempt against its task. A failing step can't take the upgrader down
            with it - only a ctrl-c gets out (recorded as a failure, then re-raised) """
        try:
            result = migration.run()
        except EOFError:  # a prompt with no stdin - nobody to answer it
            print_red("No input available")
            result = MigrationResult.skip()
        except KeyboardInterrupt:
            self._record(migration, MigrationResult.failure("Interrupted"))
            raise
        except SystemExit as e:
            if e.code in (None, 0):
                result = MigrationResult.success()
            else:
                result = MigrationResult.failure(f"Exited with code {e.code}")
        except Exception as e:
            traceback.print_exc()
            result = MigrationResult.failure(f"{type(e).__name__}: {e}")
        finally:
            # Don't let a connection a step broke leak into the next one (a test's transaction stays open)
            for connection in connections.all(initialized_only=True):
                if not connection.in_atomic_block:
                    connection.close()

        self._record(migration, result)
        if result.status == MigrationStatus.SUCCESS:
            print_green("*** task succeeded ***")
        elif result.status == MigrationStatus.FAILURE:
            print_red("*** task failed ***")
        return result

    def _record(self, migration: SubMigration, result: MigrationResult):
        if migration.task_id and result.status != MigrationStatus.SKIP:
            record_attempt(migration.task_id, success=result.status == MigrationStatus.SUCCESS,
                           note=result.note, version=self.git_version)

    def run_selection(self, migrations: list[SubMigration], then: str, stop_on_failure: bool = False) -> list[SubMigration]:
        """ Runs migrations in turn and returns those that failed. Carries on past a failure unless
            stop_on_failure, skipping any step gated after one that failed. A pull that brings new code
            restarts the process, which runs the rest of the selection and then 'then' """
        failed_ids = set()
        failed = []
        for i, migration in enumerate(migrations):
            if failed_ids and (waiting_on := self._failed_prerequisites(migration, failed_ids)):
                print_yellow(f"Skipping '{migration}' - prerequisite failed: {', '.join(waiting_on)}")
                continue
            result = self.run_step(migration)
            if result.status == MigrationStatus.FAILURE:
                failed.append(migration)
                if migration.task_id:
                    failed_ids.add(migration.task_id)
                if stop_on_failure:
                    break
            if result.new_code:
                self.restart_after_pull([m.ref for m in migrations[i + 1:]], then=then,
                                        stop_on_failure=stop_on_failure)
        return failed

    @staticmethod
    def _failed_prerequisites(migration: SubMigration, failed_ids: set[str]) -> list[str]:
        return [task_id for gate in (migration.requires or []) if gate.startswith(AFTER_PREFIX)
                and (task_id := gate.removeprefix(AFTER_PREFIX)) in failed_ids]

    def restart_after_pull(self, refs: list[str], then: str, stop_on_failure: bool):
        print_light_purple("Pulled new code - restarting the upgrader to run it")
        sys.stdout.flush()
        connections.close_all()
        resume = json.dumps({"refs": refs, "then": then, "stop_on_failure": stop_on_failure})
        manage_py = os.path.join(settings.BASE_DIR, "manage.py")
        os.execv(sys.executable, [sys.executable, manage_py, "upgrade", "--resume", resume])

    def resolve_refs(self, refs: list[str]) -> list[SubMigration]:
        """ Steps saved by restart_after_pull. Manual tasks resolve by task id, since their numbering can change """
        by_ref = {migration.ref: migration for migration in self.menu_migrations()}
        migrations = []
        for ref in refs:
            if migration := by_ref.get(ref):
                migrations.append(migration)
            else:
                print_yellow(f"'{ref}' is no longer outstanding - skipping")
        return migrations

    def resume(self, resume_json: str) -> int:
        resume = json.loads(resume_json)
        failed = self.run_selection(self.resolve_refs(resume["refs"]), then=resume["then"],
                                    stop_on_failure=resume["stop_on_failure"])
        return self.finish(resume["then"], failed)

    def finish(self, then: str, failed: list[SubMigration]) -> int:
        """ What happens after a selection: 'exit' (with a summary), 'quick' (exit if nothing is left,
            otherwise the menu) or 'menu' """
        if then == "exit":
            return self.report_failed(failed)
        if then == "quick":
            outstanding = outstanding_tasks_json()
            if not failed and not outstanding:
                print_light_purple("Quick migration was successful")
                print_restart_reminder()
                return 0
            if outstanding:
                print_red("Outstanding custom migrations, remember you can mark them all as skipped using VGs version page")
        return self.prompt()

    @staticmethod
    def report_failed(failed: list[SubMigration]) -> int:
        if failed:
            print_red("Failed steps:")
            for migration in failed:
                print_red(f"    {migration}")
            return 1
        return 0

    def quick(self) -> int:
        print_purple("-- Attempting automatic update --")
        failed = self.run_selection(self.standard_migrations(), then="quick", stop_on_failure=True)
        return self.finish("quick", failed)

    def run_steps(self, selection: str) -> int:
        failed = self.run_selection(parse_selection(selection, self.menu_migrations()), then="exit")
        return self.finish("exit", failed)

    def auto_manage(self) -> list[dict]:
        """ Auto-run every 'manage' task whose gates are satisfied, re-evaluating between passes so
            task-ordering gates (after:<task_id>) unblock as their prerequisites complete. Carries on past a
            failure; leaves blocked manage tasks and all 'other' tasks for a human. Returns the failed tasks """
        print_purple("-- Auto-running unblocked manage.py steps --")

        def fetch_runnable():
            return [task for task in outstanding_tasks_json() if task["runnable"]]

        def run_task(task) -> bool:
            return self.run_step(Upgrader.subcommand_for_json(task)).status == MigrationStatus.SUCCESS

        _, failed = run_scheduler(fetch_runnable, run_task, key=lambda task: task["id"])
        if failed:
            print_red("Failed manage.py steps (anything gated after them is still blocked):")
            for task in failed:
                print_red(f"    python3 manage.py {task['line']}")
        self.report_outstanding()
        return failed

    @staticmethod
    def report_outstanding():
        tasks = outstanding_tasks_json()
        manage = [t for t in tasks if t["category"] == "manage"]
        obsolete = [t for t in manage if not t["command_exists"]]
        manage_blocked = [t for t in manage if t["command_exists"] and not t["runnable"]]
        others = [t for t in tasks if t["category"] != "manage"]

        if manage_blocked:
            print_yellow("Blocked manage.py steps (prerequisite gate not satisfied):")
            for t in manage_blocked:
                gates = ", ".join(t["blocked_by"])
                print_yellow(f"    python3 manage.py {t['line']}    [waiting on: {gates}]")
        if obsolete:
            print_yellow("Obsolete manage steps (command no longer exists):")
            for t in obsolete:
                print_yellow(f"    {t['line']}")
        if others:
            print_yellow("Manual steps still requiring you (not auto-run):")
            for t in others:
                print_yellow(f"    {t['line']}")
        if not manage_blocked and not obsolete and not others:
            print_green("No outstanding manual steps remain.")
        if manage_blocked:
            print_purple("Satisfy a manual gate with: python3 manage.py manual_gate --satisfy <gate>")
            print_purple("List gate status with:      python3 manage.py manual_gate")
        if obsolete:
            print_purple("Mark an obsolete step complete by selecting it in the ./scripts/upgrade.sh menu")

    def prompt(self) -> int:
        while True:
            migrations = self.menu_migrations()
            self.print_menu(migrations)
            try:
                selection = input("\033[95mPlease enter a selection (keys or ranges, e.g. 'm c 2-5'): \033[00m").strip()
            except EOFError:
                print_red("No input available for the menu - exiting")
                return 1

            try:
                if selection == "q":
                    print_restart_reminder()
                    return 0
                if selection == "a":
                    self.run_selection(self.standard_migrations(), then="menu", stop_on_failure=True)
                elif selection == "am":
                    self.auto_manage()
                else:
                    try:
                        selected = parse_selection(selection, migrations)
                    except ValueError as e:
                        print_red(str(e))
                        continue
                    self.run_selection(selected, then="menu")
            except KeyboardInterrupt:
                print_red("\nInterrupted")

    @staticmethod
    def print_menu(migrations: list[SubMigration]):
        print_purple("-- Welcome to variantgrid upgrader --")
        print(f"a: automate standard steps ({Upgrader.STANDARD_KEYS}), stopping at the first failure")
        print("am: auto-run all unblocked manage.py steps (skips gated + non-manage manual steps)")
        for migration in migrations:
            if migration.key == "1":
                print("****** SPECIAL STEPS ******")
            tag = migration.status_tag()
            prefix = f"{color(YELLOW, tag)} " if tag else ""
            print(f"{migration.key}: {prefix}{migration!s}")
            if status := migration.status_line():
                print_yellow(f"    ⧗ {status}")
            for note in migration.notes or []:
                print(f"    {note}")
        print("q: exit")
