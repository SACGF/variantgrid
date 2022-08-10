from typing import Optional, List, Union, Callable

from django.db.migrations.operations.base import Operation


class ManualOperation(Operation):

    reversable = False
    reduces_to_sql = False

    def __init__(self, task_id: str, note: Optional[str] = None, test: Callable = None, *args, **kwargs):
        """ test - optional callable, only create manual operation if test returns True  """
        self.task_id = task_id
        self.note = note
        self.test = test

    def deconstruct(self):
        kwargs = {
            'task_id': self.task_id,
        }
        if self.note:
            kwargs['note'] = self.note

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
            ManualMigrationRequsted.objects.create(task=task, note=self.note)

    @staticmethod
    def operation_manage(args: List[str], note: Optional[str] = None, test: Callable = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("manage", args), note=note, test=test)

    @staticmethod
    def operation_other(args: List[str], note: Optional[str] = None, test: Callable = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("other", args), note=note, test=test)

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
    def _task_id_generate(category: str, args: List[str]) -> str:
        if isinstance(args, str):
            args = [args]
        args = [ManualOperation.escape_arg(arg) for arg in args]
        arg_string = " ".join(args)
        return f"{category}*{arg_string}"

    @staticmethod
    def task_id_manage(args: Union[str, List[str]]) -> str:
        return ManualOperation._task_id_generate(category="manage", args=args)
