"""
Branch migration_checkpoint is a dead end (#736): once its last migration is applied, every command that runs the
system checks (runserver, migrate, the upgrader's deployment_check) errors until the checkout moves to master or
the release branch, whose squashed migrations take over from here.

Entry point: check_migration_checkpoint_reached (system check).
"""
from django.conf import settings
from django.core.checks import Error, register
from django.db import DatabaseError, connection
from django.db.migrations.recorder import MigrationRecorder

CHECKPOINT_MIGRATION = ("manual", "0009_migration_checkpoint")


@register()
def check_migration_checkpoint_reached(app_configs, databases=None, **kwargs):
    if settings.UNIT_TEST:
        return []  # A test database applies the checkpoint like any other migration
    try:
        reached = MigrationRecorder(connection).migration_qs.filter(app=CHECKPOINT_MIGRATION[0],
                                                                    name=CHECKPOINT_MIGRATION[1]).exists()
    except DatabaseError:
        return []
    if not reached:
        return []
    return [Error("Checkpoint reached: switch to master / the release branch (git checkout, then scripts/upgrade.sh)",
                  id="manual.E002")]
