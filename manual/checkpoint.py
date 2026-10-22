"""
Branch migration_checkpoint is a dead end (#736): once its last migration is applied, every command that runs the
system checks (runserver, migrate, the upgrader's deployment_check) errors until the checkout moves to master or
the release branch, whose squashed migrations take over from here.

sapath moves branches with core: check_sapath_checkout_matches stops every command while the sapath checkout lacks
the migration the checkpoint depends on (sapath on its squashed master, or older than this branch). It looks at the
files rather than loading migrations, because a squashed sapath's migrations import code this branch doesn't have.

Entry points: check_migration_checkpoint_reached, check_sapath_checkout_matches (system checks).
"""
import os

from django.apps import apps
from django.conf import settings
from django.core.checks import Error, register
from django.db import DatabaseError, connection
from django.db.migrations.recorder import MigrationRecorder

CHECKPOINT_MIGRATION = ("manual", "0009_migration_checkpoint")
SAPATH_CHECKPOINT_MIGRATION = ("sapath", "0025_merge_tso500_kits_timestamps")


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


@register()
def check_sapath_checkout_matches(app_configs, **kwargs):
    if not apps.is_installed("sapath"):
        return []
    app_label, name = SAPATH_CHECKPOINT_MIGRATION
    migration_file = os.path.join(apps.get_app_config(app_label).path, "migrations", f"{name}.py")
    if os.path.exists(migration_file):
        return []
    return [Error(f"The sapath checkout has no migration {name}, which this branch's checkpoint depends on - it is "
                  "on its squashed master (or older than this branch). Check out sapath's 'migration_checkpoint' "
                  "branch too (git checkout migration_checkpoint in variantgrid_sapath).",
                  id="manual.E004")]
