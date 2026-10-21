from django.apps import AppConfig
from django.db.models.signals import pre_migrate

from manual.squashed_history import refuse_migrate_on_old_history  # also registers its system check


class ManualConfig(AppConfig):
    name = 'manual'

    def ready(self):
        pre_migrate.connect(refuse_migrate_on_old_history, sender=self)
