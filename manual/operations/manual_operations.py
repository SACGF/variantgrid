from typing import Optional, List, Union

from django.db.migrations.operations.base import Operation


class ManualOperation(Operation):

    reversable = False
    reduces_to_sql = False

    def __init__(self, task_id: str, note: Optional[str] = None, *args, **kwargs):
        self.task_id = task_id
        self.note = note

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

    def run(self, apps):
        """
        Call this to run the operation right away in RunPython
        :param apps: where we can call get_model
        """
        ManualMigrationTask = apps.get_model('manual', 'ManualMigrationTask')
        ManualMigrationRequsted = apps.get_model('manual', 'ManualMigrationRequired')

        task, _ = ManualMigrationTask.objects.get_or_create(pk=self.task_id)
        ManualMigrationRequsted.objects.create(task=task, note=self.note)

    @staticmethod
    def operation_manage(args: List[str], note:Optional[str] = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("manage", args), note=note)

    @staticmethod
    def operation_other(args: List[str], note: Optional[str] = None):
        return ManualOperation(task_id=ManualOperation._task_id_generate("other", args), note=note)

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
