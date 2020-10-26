from django.db.migrations.operations.base import Operation


class ManualOperation(Operation):

    reversable = False
    reduces_to_sql = False

    def __init__(self, app: str, release_note: str, *args, **kwargs):
        self.app = app
        self.release_note = release_note

    def deconstruct(self):
        kwargs = {
            'app': self.app,
            'release_note': self.release_note
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
        Deployment = apps.get_model('manual', 'Deployment')
        Deployment.objects.create(release_note = self.release_note, app = self.app)