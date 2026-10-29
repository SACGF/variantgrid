"""
Migration history from before the squash (#736): every app's migrations were replaced by fresh initials plus a
`seed_data` migration, with no `replaces`, so a database must reach this code through branch
`migration_checkpoint`. That branch's last migration (CHECKPOINT_MIGRATION) depends on every app's final
pre-squash migration and refuses to apply while manual tasks are outstanding, so its row in django_migrations
means the database holds the whole old history and every one-off step it registered.

A database is in one of three states, judged on the django_migrations rows of our own apps (not site-packages)
that name no migration on disk. Once NEW_HISTORY_MIGRATION is recorded (a fresh build, or the reset) only rows applied
before it are old history: one applied after it is a post-squash migration renamed since, which Django ignores.
- current: no old-history rows - a fresh database, or one already reset
- needs_checkpoint: rows from the old history and no checkpoint row - upgrade through the checkpoint branch first
- needs_reset: rows from the old history and the checkpoint row - delete them and record the new initial
  migrations as applied (every one with `initial = True`, which the squash's seed migrations also set)

Entry points: history_state, reset_squashed_history (`manage.py reset_squashed_migration_history`, run by hand
once per deployment so the rewrite of django_migrations is seen), check_squashed_history (a database system check, so
migrate, `check --database default` and deployment_check run it), refuse_migrate_on_old_history (pre_migrate) and
refuse_to_serve_old_history (gunicorn's wsgi.py and celery's worker_init), so a server left on master that pulls the
squash without going through the checkpoint stops rather than serving the new code on the old schema.
A core checkout whose sapath checkout is still pre-squash can't even build the migration graph: each of those entry
points reports CheckoutBeforeSquashError (manual.E003 as a check) naming the stale checkout.
"""
import os
from enum import Enum

from celery.signals import worker_init
from django.apps import apps
from django.core.checks import Error, Tags, register
from django.db import DatabaseError, connections, transaction
from django.db.migrations.exceptions import NodeNotFoundError
from django.db.migrations.loader import MigrationLoader
from django.db.migrations.recorder import MigrationRecorder

CHECKPOINT_BRANCH = "migration_checkpoint"
CHECKPOINT_MIGRATION = ("manual", "0009_migration_checkpoint")
NEW_HISTORY_MIGRATION = ("snpdb", "0002_seed_data")  # The old history had no migration of this name


class HistoryState(Enum):
    CURRENT = "current"
    NEEDS_CHECKPOINT = "needs_checkpoint"
    NEEDS_RESET = "needs_reset"


def _local_app_labels() -> set[str]:
    return {app_config.label for app_config in apps.get_app_configs()
            if "site-packages" not in os.path.realpath(app_config.path)}


class CheckoutBeforeSquashError(RuntimeError):
    """ Another checkout (sapath) still has its pre-squash migrations, which name ones the squash removed, so the
        migration graph can't be built at all - before any database state can be judged """


def _migration_loader(using: str) -> MigrationLoader:
    try:
        return MigrationLoader(connections[using], ignore_no_migrations=True)
    except NodeNotFoundError as e:
        missing_app = e.node[0]
        if missing_app not in _local_app_labels():
            raise
        origin_app = getattr(e.origin, "app_label", e.origin)
        raise CheckoutBeforeSquashError(
            f"'{origin_app}' migrations depend on {e.node}, which the migration squash (#736) removed: the checkout "
            f"holding '{origin_app}' is from before the squash. Update it too (core and sapath move together), "
            "then rerun.") from e


def _stale_migrations(loader: MigrationLoader) -> set[tuple[str, str]]:
    local_apps = _local_app_labels()
    stale = {key: row for key, row in loader.applied_migrations.items()
             if key[0] in local_apps and key not in loader.disk_migrations and key != CHECKPOINT_MIGRATION}
    if new_history_row := loader.applied_migrations.get(NEW_HISTORY_MIGRATION):
        return {key for key, row in stale.items() if row.applied < new_history_row.applied}
    return set(stale)


def history_state(using: str = "default") -> HistoryState:
    loader = _migration_loader(using)
    if not _stale_migrations(loader):
        return HistoryState.CURRENT
    if CHECKPOINT_MIGRATION in loader.applied_migrations:
        return HistoryState.NEEDS_RESET
    return HistoryState.NEEDS_CHECKPOINT


def history_state_message(state: HistoryState) -> str:
    if state == HistoryState.NEEDS_CHECKPOINT:
        return (f"This database has migrations from before the squash (#736). Check out branch '{CHECKPOINT_BRANCH}', "
                "run scripts/upgrade.sh there (migrate, then every outstanding manual task), then come back to this "
                "branch and run 'python3 manage.py reset_squashed_migration_history'.")
    if state == HistoryState.NEEDS_RESET:
        return ("This database has passed the migration checkpoint but still holds the old migration history: run "
                "'python3 manage.py reset_squashed_migration_history', then scripts/upgrade.sh.")
    return "Migration history is current."


def reset_squashed_history(using: str = "default") -> int:
    """ Replaces our apps' django_migrations rows with the new initial migrations. Returns the number deleted """
    if (state := history_state(using)) != HistoryState.NEEDS_RESET:
        raise ValueError(f"Only a database that has passed the checkpoint can be reset: {history_state_message(state)}")
    connection = connections[using]
    loader = _migration_loader(using)
    local_apps = _local_app_labels()
    initial_migrations = [key for key, migration in loader.disk_migrations.items()
                          if key[0] in local_apps and migration.initial]
    recorder = MigrationRecorder(connection)
    with transaction.atomic(using=using):
        old_rows = recorder.migration_qs.filter(app__in=local_apps).exclude(app=CHECKPOINT_MIGRATION[0],
                                                                            name=CHECKPOINT_MIGRATION[1])
        deleted, _ = old_rows.delete()
        for app_label, name in initial_migrations:
            recorder.record_applied(app_label, name)
    return deleted


@register(Tags.database)
def check_squashed_history(app_configs, databases=None, **kwargs):
    try:
        state = history_state()
    except DatabaseError:
        return []  # No database (a build step, or before it is created) - migrate will say so
    except CheckoutBeforeSquashError as e:
        return [Error(str(e), id="manual.E003")]
    if state == HistoryState.CURRENT:
        return []
    return [Error(history_state_message(state), id="manual.E001")]


def refuse_migrate_on_old_history(sender, using="default", **kwargs):
    """ pre_migrate handler: migrating an old-history database would treat a new 0001_initial as applied because
        it shares the old one's name, then apply the rest of the new history on top of the old schema """
    state = history_state(using)
    if state != HistoryState.CURRENT:
        raise RuntimeError(history_state_message(state))


def refuse_to_serve_old_history(using: str = "default"):
    try:
        state = history_state(using)
    except DatabaseError:
        return  # Database down - the services fail on their own, as before
    finally:
        connections[using].close()  # Opened in the parent: gunicorn and celery fork workers from it
    if state != HistoryState.CURRENT:
        raise RuntimeError(history_state_message(state))


@worker_init.connect
def refuse_celery_worker_on_old_history(**kwargs):
    try:
        refuse_to_serve_old_history()
    except RuntimeError as e:
        raise SystemExit(str(e)) from e  # celery logs and carries on past an Exception from a signal handler
