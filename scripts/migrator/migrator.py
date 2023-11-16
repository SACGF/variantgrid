import json
import os
import socket
import subprocess
import sys
import traceback
from enum import Enum, auto
from functools import partial
from subprocess import Popen
from typing import List, Optional, Dict, Callable

import requests

from library.git import Git

COMMAND_ALIASES = {
    "python": "python3"
}


def print_color(color: str, skk: str):
    print((color + "{}\033[00m").format(skk))


print_red = partial(print_color, "\033[91m")
print_green = partial(print_color, "\033[92m")
print_yellow = partial(print_color, "\033[93m")
print_light_purple = partial(print_color, "\033[94m")
print_purple = partial(print_color, "\033[95m")
print_cyan = partial(print_color, "\033[96m")
print_light_gray = partial(print_color, "\033[97m")
print_black = partial(print_color, "\033[98m")


def substitute_aliases(args: List[str]):
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

    def using(self, key: Optional[str] = None, task_id: Optional[str] = None, notes: Optional[List[str]] = None):
        if key:
            self.key = key
        if task_id:
            self.task_id = task_id
        if notes:
            self.notes = notes
        return self

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
            selection = input("\033[95m{}\033[00m".format("Please enter a selection: "))
            selection = selection.strip().lower()
            if selection == "y":
                return MigrationResult.success()
            if selection == "n":
                return MigrationResult.failure()
            if selection == "x":
                return MigrationResult.skip()
            print(f"Unexpected input - \"{selection}\"")


class GitSubMigration(SubMigration):

    def __init__(self, args: List[str]):
        super().__init__()
        self.args = args

    @property
    def effective_args(self):
        return substitute_aliases(self.args)

    def __str__(self):
        return "git " + " ".join(self.effective_args)

    def run(self) -> MigrationResult:
        print_cyan(str(self))
        completed_process = subprocess.run("git pull", shell=True, check=False)
        if completed_process.returncode != 0:
            return MigrationResult(status=MigrationStatus.FAILURE,
                                   note=f"Subprocess failed with error code {completed_process.returncode}")
        return MigrationResult(status=MigrationStatus.SUCCESS)


class CommandSubMigration(SubMigration):

    def __init__(self, args: List[str]):
        super().__init__()
        self.args = args

    @property
    def effective_args(self) -> List[str]:
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
    def manage_py(args: List[str]):
        use_args = ["python", "manage.py"]
        use_args.extend(args)
        return CommandSubMigration(args=use_args)

    @staticmethod
    def bash(args: List[str]):
        return CommandSubMigration(args=args)

    @staticmethod
    def git(args: List[str]):
        return GitSubMigration(args=args)

    @staticmethod
    def script(args: List[str]):
        args[0] = f"./scripts/{args[0]}"
        return CommandSubMigration(args)

    @staticmethod
    def python(the_method: Callable):
        return PythonSubMigration(the_method)


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
            response = requests.post("https://api.rollbar.com/api/1/deploy/", data=data)

            if response.status_code == 200:
                print("Deployment recorded in Rollbar.")
                subprocess.run(["python", "manage.py", "deployed"])
            else:
                print(f"Failed to record deployment in Rollbar. Response: {response.text}")
        except Exception as e:
            print(f"Error recording deployment in Rollbar: {str(e)}")

        return MigrationResult.success()

    STANDARD_MIGRATIONS: List[SubMigration] = [
        CommandSubMigration.git(["pull"]).using(key="g", task_id="git*pull"),
        CommandSubMigration.manage_py(["migrate"]).using(key="m", task_id="manage*migrate"),
        CommandSubMigration.manage_py(["collectstatic_js_reverse"]).using(key="r",
                                                                          task_id="manage*collectstatic_js_reverse"),
        # collectstatic without warning for conflicting files has been an issue for 6 years
        # see https://code.djangoproject.com/ticket/26583 maybe it'll get fixed soon? For now (since we've never had
        # a problem) just turn off all verbosity
        CommandSubMigration.manage_py(["collectstatic", "-v", "0", "--noinput"]).using(key="c",
                                                                                       task_id="manage*collectstatic"),
    ]

    @property
    def standard_migrations(self):
        return Migrator.STANDARD_MIGRATIONS + [CommandSubMigration.python(self.notify_deployed).using(key="d")]

    def __init__(self):
        self.migrations = self.standard_migrations
        self.has_custom_migrations = False
        self.git_version = Migrator.get_git_ver()
        self.rollbar_token = None

    def refresh_migrations(self):
        try:
            migrations: List[SubMigration] = self.standard_migrations

            command = substitute_aliases(["python", "manage.py", "manual_outstanding"])
            self.has_custom_migrations = False
            with Popen(command, stdout=subprocess.PIPE) as proc:
                # stdout_text = proc.stdout.read()
                task_json = json.load(proc.stdout)
                if token := task_json.get('ROLLBAR_ACCESS_TOKEN'):
                    self.rollbar_token = token
                command = 1
                for task in task_json["tasks"]:
                    self.has_custom_migrations = True
                    migrations.append(Migrator.subcommand_for_json(task).using(key=f"{command}"))
                    command += 1

            self.migrations = migrations
        except:
            print("Unable to retrieve outstanding commands")
            traceback.print_exc()

    def record_attempt(self, task_id: str, success: bool = True):
        args = ["python", "manage.py", "manual_complete", "--id", task_id, "--ver", self.git_version]
        if not success:
            args.extend(["--failed"])
        command = substitute_aliases(args)
        with Popen(command, stdout=subprocess.PIPE) as proc:
            proc.stdout.read()  # do we care about this output

    @staticmethod
    def get_git_ver() -> str:
        return Git().hash

    @staticmethod
    def subcommand_for_json(task: Dict) -> SubMigration:
        task_id = task["id"]
        category = task["category"]
        line = task["line"]
        notes = task.get("notes")
        if category == "manage":
            args = ["python", "manage.py", line]
            return CommandSubMigration(args).using(task_id=task_id, notes=notes)
        return ManualSubMigration(line).using(task_id=task_id, notes=notes)

    def prompt(self, refresh: bool = True):
        if refresh:
            self.refresh_migrations()
        keys = []
        print_purple("-- Welcome to variantgrid upgrader --")
        print("a: automate standard steps (runs git, migrate, collectstatic_js_reverse, collectstatic, deployed)")
        for migration in self.migrations:
            if migration.key == "1":
                print("****** SPECIAL STEPS ******")
            print(f"{migration.key}: {str(migration)}")
            keys.append(migration.key)
            if migration.notes:
                for note in migration.notes:
                    print(f"    {note}")
        print("q: exit")
        keys.append("q")
        keys.append("a")

        selected_migration: Optional[SubMigration] = None
        while selected_migration is None:
            selection = input("\033[95m{}\033[00m".format(f"Please enter a selection ({','.join(keys)}): "))
            selection = selection.strip()

            if selection == "a":
                self.run_and_re_prompt(self.standard_migrations)
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

    def run_and_re_prompt(self, migrations: List[SubMigration]):
        self.run_and_callback(migrations, lambda _: self.prompt())

    def run_and_callback(self, migrations: List[SubMigration], callback: Callable[[bool], None]):
        migration: Optional[SubMigration] = None
        if migrations:
            migration = migrations.pop(0)

        if migration is None:
            callback(True)
        else:
            result = migration.run()
            if result.status != MigrationStatus.SKIP:
                if migration.task_id:
                    self.record_attempt(task_id=migration.task_id, success=result.status == MigrationStatus.SUCCESS)
            if result.status == MigrationStatus.SUCCESS:
                print_green("*** task succeeded ***")
                self.run_and_callback(migrations, callback)
            else:
                # don't continue to auto-run migration step on failure
                print_red("*** task failed ***")
                callback(False)


if __name__ == '__main__':
    quick_mode = False
    if len(sys.argv) > 1:
        quick_mode = sys.argv[1] == '--quick'

    migrator = Migrator()
    if quick_mode:
        print_purple("-- Attempting automatic update --")
        migrator.run_and_quit_if_success()
    else:
        migrator.prompt()
