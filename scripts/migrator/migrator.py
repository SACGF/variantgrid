import json
import os
import socket
import subprocess
import sys
import traceback
from collections.abc import Callable
from enum import Enum, auto
from functools import partial
from subprocess import Popen
from typing import Optional

import requests

from library.constants import MINUTE_SECS
from library.git import Git

COMMAND_ALIASES = {
    "python": "python3"
}


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
print_light_gray = partial(print_color, "\033[97m")
print_black = partial(print_color, "\033[98m")


def substitute_aliases(args: list[str]):
    return [COMMAND_ALIASES.get(arg, arg) for arg in args]


class MigrationStatus(Enum):
    SUCCESS = auto()
    FAILURE = auto()
    SKIP = auto()


class MigrationResult:

    def __init__(self, status: MigrationStatus, note: Optional[str] = None):
        self.status = status
        self.note = note

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
        self.blocked_by = None       # gate names not yet satisfied (see manual.gates)
        self.command_exists = None   # False for an obsolete 'manage' task whose command was deleted

    def using(self, key: Optional[str] = None, task_id: Optional[str] = None, notes: Optional[list[str]] = None):
        if key:
            self.key = key
        if task_id:
            self.task_id = task_id
        if notes:
            self.notes = notes
        return self

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

    def __init__(self, the_method: Callable[[], MigrationResult]):
        super().__init__()
        self.the_method = the_method

    def run(self):
        return self.the_method()

    def __str__(self):
        return self.the_method.__name__


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


class GitSubMigration(SubMigration):
    """ Pulls, then installs any dependency changes that came down with it. These are one step as
        the migrator calls manage.py between steps (which imports every 3rd party package), so the
        new requirements have to be in place before we hand back. """

    def __init__(self, args: list[str]):
        super().__init__()
        self.args = args

    @property
    def effective_args(self):
        return substitute_aliases(self.args)

    def __str__(self):
        return "git " + " ".join(self.effective_args) + " + install requirements"

    def run(self) -> MigrationResult:
        print_cyan(str(self))
        for command in ["git pull", "./scripts/install_requirements.sh"]:
            completed_process = subprocess.run(command, shell=True, check=False)
            if completed_process.returncode != 0:
                return MigrationResult(status=MigrationStatus.FAILURE,
                                       note=f"'{command}' failed with error code {completed_process.returncode}")
        return MigrationResult(status=MigrationStatus.SUCCESS)


class CommandSubMigration(SubMigration):

    def __init__(self, args: list[str]):
        super().__init__()
        self.args = args

    @property
    def effective_args(self) -> list[str]:
        return substitute_aliases(self.args)

    def run(self) -> MigrationResult:
        print_cyan(str(self))
        print_purple("-----------")
        cmd = str(self)
        completed_process = subprocess.run(cmd, shell=True, check=False)
        print_purple("-----------")
        if completed_process.returncode != 0:
            note = f"Subprocess failed with error code {completed_process.returncode}"
            return MigrationResult(status=MigrationStatus.FAILURE, note=note)
        return MigrationResult(status=MigrationStatus.SUCCESS)

    def __str__(self):
        return " ".join(self.effective_args)

    @staticmethod
    def manage_py(args: list[str]):
        use_args = ["python", "manage.py"]
        use_args.extend(args)
        return CommandSubMigration(args=use_args)

    @staticmethod
    def bash(args: list[str]):
        return CommandSubMigration(args=args)

    @staticmethod
    def git(args: list[str]):
        return GitSubMigration(args=args)

    @staticmethod
    def script(args: list[str]):
        args[0] = f"./scripts/{args[0]}"
        return CommandSubMigration(args)

    @staticmethod
    def python(the_method: Callable):
        return PythonSubMigration(the_method)


def run_scheduler(fetch_runnable, run_task, max_passes: int = 1000):
    """ Pure scheduler loop, decoupled from Django/subprocess for testability.

        Repeatedly asks fetch_runnable() for the currently-runnable tasks and runs each, re-fetching
        after every pass so ordering gates (after:<task_id>) release as their prerequisites complete.
        Stops on the first failure. max_passes is a runaway guard.

        fetch_runnable() -> list[task];  run_task(task) -> bool (True on success).
        Returns (ran, failed): tasks run in order, and the task that failed (or None). """
    ran = []
    for _ in range(max_passes):
        runnable = fetch_runnable()
        if not runnable:
            return ran, None
        for task in runnable:
            ran.append(task)
            if not run_task(task):
                return ran, task
    return ran, None


class Migrator:

    def notify_deployed(self):
        if not self.rollbar_token:
            self.refresh_migrations()
            if not self.rollbar_token:
                print_red("No rollbar token found")
                return MigrationResult.failure()

        environment = socket.gethostname().lower().split('.')[0].replace('-', '')
        local_username = os.getenv("USER", "")
        revision = subprocess.check_output(["git", "rev-parse", "--verify", "HEAD"]).decode("utf-8").strip()

        data = {
            "access_token": self.rollbar_token,
            "environment": environment,
            "revision": revision,
            "local_username": local_username,
        }

        try:
            response = requests.post("https://api.rollbar.com/api/1/deploy/", data=data, timeout=MINUTE_SECS)

            if response.status_code == 200:
                subprocess.run(substitute_aliases(["python", "manage.py", "deployed"]))
            else:
                print(f"Failed to record deployment in Rollbar. Response: {response.text}")
        except Exception as e:
            print(f"Error recording deployment in Rollbar: {e!s}")

        return MigrationResult.success()

    STANDARD_MIGRATIONS: list[SubMigration] = [
        CommandSubMigration.git(["pull"]).using(key="g", task_id="git*pull"),
        CommandSubMigration.manage_py(["migrate"]).using(key="m", task_id="manage*migrate"),
        CommandSubMigration.manage_py(["collectstatic_js_reverse"]).using(key="r",
                                                                          task_id="manage*collectstatic_js_reverse"),
        # collectstatic without warning for conflicting files has been an issue for 6 years
        # see https://code.djangoproject.com/ticket/26583 maybe it'll get fixed soon? For now (since we've never had
        # a problem) just turn off all verbosity
        # --clear so a moved file whose mtime looks unchanged doesn't leave a stale copy behind
        CommandSubMigration.manage_py(["collectstatic", "-v", "0", "--noinput", "--clear"]).using(key="c",
                                                                                       task_id="manage*collectstatic"),
        CommandSubMigration.manage_py(["deployment_check", "--die-if-invalid", "--quiet"]).using(key="k",
                                                                                       task_id="manage*deployment_check"),
    ]

    @property
    def standard_migrations(self):
        return Migrator.STANDARD_MIGRATIONS + [CommandSubMigration.python(self.notify_deployed).using(key="d")]

    def __init__(self):
        self.migrations = self.standard_migrations
        self.has_custom_migrations = False
        self.git_version = Migrator.get_git_ver()
        self.rollbar_token = None

    def _fetch_outstanding(self) -> dict:
        """ Call the manual_outstanding management command and return its parsed JSON envelope.
            Each task carries 'category', 'line', 'requires', 'blocked_by' and 'runnable'. """
        command = substitute_aliases(["python", "manage.py", "manual_outstanding"])
        with Popen(command, stdout=subprocess.PIPE) as proc:
            task_json = json.load(proc.stdout)
        if token := task_json.get('ROLLBAR_ACCESS_TOKEN'):
            self.rollbar_token = token
        return task_json

    def refresh_migrations(self):
        try:
            migrations: list[SubMigration] = self.standard_migrations

            task_json = self._fetch_outstanding()
            self.has_custom_migrations = False
            command = 1
            for task in task_json["tasks"]:
                self.has_custom_migrations = True
                migrations.append(Migrator.subcommand_for_json(task).using(key=f"{command}"))
                command += 1

            self.migrations = migrations
        except Exception:
            print("Unable to retrieve outstanding commands")
            traceback.print_exc()

    def run_auto_manage(self):
        """ Auto-run every 'manage' task whose gates are satisfied, re-evaluating between passes so
            task-ordering gates (after:<task_id>) unblock as their prerequisites complete. Stops on
            the first failure; leaves blocked manage tasks and all 'other' tasks for a human. """
        print_purple("-- Auto-running unblocked manage.py steps --")

        def fetch_runnable():
            try:
                task_json = self._fetch_outstanding()
            except Exception:
                print_red("Unable to retrieve outstanding commands")
                traceback.print_exc()
                return []
            return [task for task in task_json.get("tasks", []) if task.get("runnable")]

        def run_task(task) -> bool:
            migration = Migrator.subcommand_for_json(task)
            result = migration.run()
            success = result.status == MigrationStatus.SUCCESS
            if migration.task_id:
                self.record_attempt(task_id=migration.task_id, success=success, note=result.note)
            if success:
                print_green("*** task succeeded ***")
            else:
                print_red("*** task failed - stopping auto-run ***")
            return success

        run_scheduler(fetch_runnable, run_task)
        self.report_outstanding()

    def report_outstanding(self):
        try:
            task_json = self._fetch_outstanding()
        except Exception:
            return
        tasks = task_json.get("tasks", [])
        manage = [t for t in tasks if t["category"] == "manage"]
        obsolete = [t for t in manage if not t.get("command_exists")]
        manage_blocked = [t for t in manage if t.get("command_exists") and not t.get("runnable")]
        others = [t for t in tasks if t["category"] != "manage"]

        if manage_blocked:
            print_yellow("Blocked manage.py steps (prerequisite gate not satisfied):")
            for t in manage_blocked:
                gates = ", ".join(t.get("blocked_by") or [])
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

    def record_attempt(self, task_id: str, success: bool = True, note: Optional[str] = None):
        args = ["python", "manage.py", "manual_complete", "--id", task_id, "--ver", self.git_version]
        if not success:
            args.extend(["--failed"])
        if note:
            args.extend(["--note", note])
        command = substitute_aliases(args)
        with Popen(command, stdout=subprocess.PIPE) as proc:
            proc.stdout.read()  # do we care about this output

    @staticmethod
    def get_git_ver() -> str:
        return Git().hash

    @staticmethod
    def subcommand_for_json(task: dict) -> SubMigration:
        task_id = task["id"]
        category = task["category"]
        line = task["line"]
        notes = task.get("notes")
        if category == "manage":
            if task.get("command_exists") is False:
                sub_migration = ObsoleteSubMigration(line).using(task_id=task_id, notes=notes)
            else:
                args = ["python", "manage.py", line]
                sub_migration = CommandSubMigration(args).using(task_id=task_id, notes=notes)
        else:
            sub_migration = ManualSubMigration(line).using(task_id=task_id, notes=notes)
        sub_migration.blocked_by = task.get("blocked_by")
        sub_migration.command_exists = task.get("command_exists")
        return sub_migration

    def prompt(self, refresh: bool = True):
        if refresh:
            self.refresh_migrations()
        keys = []
        print_purple("-- Welcome to variantgrid upgrader --")
        print("a: automate standard steps (runs git + install requirements, migrate, collectstatic_js_reverse, collectstatic, deployment_check, deployed)")
        print("am: auto-run all unblocked manage.py steps (skips gated + non-manage manual steps)")
        for migration in self.migrations:
            if migration.key == "1":
                print("****** SPECIAL STEPS ******")
            tag = migration.status_tag()
            prefix = f"{color(YELLOW, tag)} " if tag else ""
            print(f"{migration.key}: {prefix}{migration!s}")
            keys.append(migration.key)
            if status := migration.status_line():
                print_yellow(f"    ⧗ {status}")
            if migration.notes:
                for note in migration.notes:
                    print(f"    {note}")
        print("q: exit")
        keys.append("q")
        keys.append("a")
        keys.append("am")

        selected_migration: Optional[SubMigration] = None
        while selected_migration is None:
            selection = input(f"\033[95mPlease enter a selection ({','.join(keys)}): \033[00m")
            selection = selection.strip()

            if selection == "a":
                self.run_and_re_prompt(self.standard_migrations)
                return

            if selection == "am":
                self.run_auto_manage()
                self.prompt()
                return

            if selection == "q":
                sys.exit(0)

            selected_options = [migration for migration in self.migrations if migration.key == selection]
            if selected_options:
                selected_migration = selected_options[0]

        self.run_and_re_prompt([selected_migration])

    def run_and_quit_if_success(self):
        def on_complete(success: bool):
            self.refresh_migrations()
            if success:
                if not self.has_custom_migrations:
                    print_light_purple("Quick migration was successful")
                    sys.exit(0)
            if self.has_custom_migrations:
                print_red("Outstanding custom migrations, remember you can mark them all as skipped using VGs version page")

            self.prompt(refresh=False)

        self.run_and_callback(migrations=self.standard_migrations, callback=on_complete)

    def run_and_re_prompt(self, migrations: list[SubMigration]):
        self.run_and_callback(migrations, lambda _: self.prompt())

    def run_and_callback(self, migrations: list[SubMigration], callback: Callable[[bool], None]):
        migration: Optional[SubMigration] = None
        if migrations:
            migration = migrations.pop(0)

        if migration is None:
            callback(True)
        else:
            result = migration.run()
            if result.status != MigrationStatus.SKIP:
                if migration.task_id:
                    self.record_attempt(task_id=migration.task_id,
                                        success=result.status == MigrationStatus.SUCCESS,
                                        note=result.note)
            if result.status == MigrationStatus.SUCCESS:
                print_green("*** task succeeded ***")
                self.run_and_callback(migrations, callback)
            else:
                # don't continue to auto-run migration step on failure
                print_red("*** task failed ***")
                callback(False)


USAGE = """variantgrid upgrader

Usage: ./scripts/upgrade.sh [MODE]
       python3 scripts/migrator/migrator.py [MODE]

Modes:
  (no argument)    Interactive menu (also offers 'a' and 'am' below).
  --quick          Run the standard steps and quit if nothing else is outstanding:
                   git pull + install requirements, migrate, collectstatic_js_reverse,
                   collectstatic, deployment_check, deployed.
  --auto-manage    Plough through all outstanding manual steps automatically: run every
                   unblocked 'manage' step (re-evaluating between passes so ordering gates
                   release as prerequisites finish), stopping on the first failure. Gated,
                   obsolete, and non-manage human steps are reported, not run. Select an
                   obsolete step in the interactive menu to mark it complete.
  --help, -h       Show this help.

Environment:
  VG_INSTALL_REQUIREMENTS=0    Skip installing requirements, so an install that won't work
                               (no network, a package that won't build) can't hold up the
                               rest of the upgrade. Dependency changes in this upgrade will
                               be missing, so migrate/collectstatic may fail.
"""


if __name__ == '__main__':
    arg = sys.argv[1] if len(sys.argv) > 1 else None

    if arg in ('--help', '-h'):
        print(USAGE)
        sys.exit(0)

    migrator = Migrator()
    if arg == '--quick':
        print_purple("-- Attempting automatic update --")
        migrator.run_and_quit_if_success()
    elif arg == '--auto-manage':
        migrator.run_auto_manage()
    else:
        migrator.prompt()
