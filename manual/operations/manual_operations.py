from collections.abc import Callable
from typing import Optional, Union

from django.db.migrations.operations.base import Operation


class ManualOperation(Operation):

    reversable = False
    reduces_to_sql = False

    def __init__(self, task_id: str, note: Optional[str] = None, test: Callable = None,
                 requires: Optional[list[str]] = None):
        """ test - optional callable, only create manual operation if test returns True
            requires - gate names (see manual.gates) that must be satisfied before this task
                       may be auto-run by the migrator """
        self.task_id = task_id
        self.note = note
        self.test = test
        self.requires = requires

    def deconstruct(self):
        kwargs = {
            'task_id': self.task_id,
        }
        if self.note:
            kwargs['note'] = self.note
        if self.requires:
            kwargs['requires'] = self.requires

        return (
            self.__class__.__qualname__,
            [],
            kwargs
        )

    def describe(self):
        return f"Registering manual migration : {self.task_id}"

    def state_forwards(self, app_label, state):
        pass

    def database_forwards(self, app_label, schema_editor, from_state, to_state):
        self.run(to_state.apps)

    def database_backwards(self, app_label, schema_editor, from_state, to_state):
        self.run(to_state.apps, reverse=True)

    def run(self, apps, reverse=False):
        """
        Call this to run the operation right away in RunPython
        :param apps: where we can call get_model
        :param reverse: for when migrating backwards (clean up)
        """
        if self.test:
            if not self.test(apps):
                return  # Skip

        ManualMigrationTask = apps.get_model('manual', 'ManualMigrationTask')
        ManualMigrationRequsted = apps.get_model('manual', 'ManualMigrationRequired')

        if reverse:
            ManualMigrationTask.objects.filter(pk=self.task_id).delete()
        else:
            task, _ = ManualMigrationTask.objects.get_or_create(pk=self.task_id)
            # 'requires' persists only once the field exists (historical migration states that
            # predate it fall back to the one-off backfill migration).
            if self.requires and hasattr(task, "requires"):
                task.requires = list(self.requires)
                task.save(update_fields=["requires"])
            if callable(self.note):
                note = self.note(apps)
            else:
                note = self.note
            ManualMigrationRequsted.objects.create(task=task, note=note)

    @staticmethod
    def operation_manage(args: list[str], note: Optional[str] = None, test: Callable = None,
                         requires: Optional[list[str]] = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("manage", args), note=note, test=test,
                               requires=requires)

    @staticmethod
    def operation_other(args: list[str], note: Optional[str] = None, test: Callable = None,
                        requires: Optional[list[str]] = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("other", args), note=note, test=test,
                               requires=requires)

    @staticmethod
    def escape_arg(arg: str):
        quote = False
        if "\"" in arg:
            quote = True
            arg = arg.replace("\"", "\\\"")
        if " " in arg:
            quote = True
        if quote:
            arg = f"\"{arg}\""
        return arg

    @staticmethod
    def _task_id_generate(category: str, args: list[str]) -> str:
        if isinstance(args, str):
            args = [args]
        args = [ManualOperation.escape_arg(arg) for arg in args]
        arg_string = " ".join(args)
        return f"{category}*{arg_string}"

    @staticmethod
    def task_id_manage(args: Union[str, list[str]]) -> str:
        return ManualOperation._task_id_generate(category="manage", args=args)

    @staticmethod
    def task_id_other(args: Union[str, list[str]]) -> str:
        return ManualOperation._task_id_generate(category="other", args=args)
