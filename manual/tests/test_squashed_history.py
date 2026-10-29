from datetime import timedelta
from unittest.mock import patch

from django.core.management import CommandError, call_command
from django.db import connection
from django.db.migrations.exceptions import NodeNotFoundError
from django.db.migrations.recorder import MigrationRecorder
from django.test import TestCase

from manual.squashed_history import (
    CHECKPOINT_MIGRATION,
    NEW_HISTORY_MIGRATION,
    CheckoutBeforeSquashError,
    HistoryState,
    history_state,
    refuse_celery_worker_on_old_history,
    reset_squashed_history,
)

OLD_MIGRATION = ("snpdb", "0297_retire_somalier_download_manual_task")


class SquashedHistoryTest(TestCase):

    def setUp(self):
        self.recorder = MigrationRecorder(connection)

    def _record_old_history(self):
        self.recorder.migration_qs.filter(app="snpdb", name__in=["0001_initial", "0002_seed_data"]).delete()
        self.recorder.record_applied(*OLD_MIGRATION)

    def test_fresh_database_is_current(self):
        self.assertEqual(history_state(), HistoryState.CURRENT)

    def test_renamed_post_squash_migration_is_current(self):
        """ A fresh build can apply a migration before NEW_HISTORY_MIGRATION; renamed since, its row is not old history """
        new_history_row = self.recorder.migration_qs.get(app=NEW_HISTORY_MIGRATION[0], name=NEW_HISTORY_MIGRATION[1])
        self.recorder.migration_qs.create(app="eventlog", name="0002_renamed_since",
                                          applied=new_history_row.applied - timedelta(seconds=1))
        self.assertEqual(history_state(), HistoryState.CURRENT)

    def test_old_history_without_checkpoint_refuses_reset(self):
        self._record_old_history()
        self.assertEqual(history_state(), HistoryState.NEEDS_CHECKPOINT)
        with self.assertRaises(CommandError):
            call_command("reset_squashed_migration_history")

    def test_reset_swaps_old_rows_for_initials(self):
        self._record_old_history()
        self.recorder.record_applied(*CHECKPOINT_MIGRATION)
        self.assertEqual(history_state(), HistoryState.NEEDS_RESET)

        reset_squashed_history()

        self.assertEqual(history_state(), HistoryState.CURRENT)
        applied = self.recorder.applied_migrations()
        self.assertNotIn(OLD_MIGRATION, applied)
        self.assertIn(("snpdb", "0001_initial"), applied)
        self.assertIn(NEW_HISTORY_MIGRATION, applied)
        self.assertIn(CHECKPOINT_MIGRATION, applied)

    def test_celery_worker_exits_on_old_history(self):
        """ celery swallows an Exception from a worker_init handler, so only SystemExit stops the worker """
        with patch("manual.squashed_history.refuse_to_serve_old_history", side_effect=RuntimeError("old history")):
            with self.assertRaises(SystemExit):
                refuse_celery_worker_on_old_history()

    def test_checkout_from_before_squash_names_it(self):
        """ Only a missing migration of one of our apps means another checkout predates the squash """
        def node_not_found(node):
            return NodeNotFoundError("missing", node, origin=type("Migration", (), {"app_label": "sapath"})())

        with patch("manual.squashed_history.MigrationLoader", side_effect=node_not_found(("annotation", "0002_old"))):
            with self.assertRaisesMessage(CheckoutBeforeSquashError, "'sapath' migrations depend on"):
                history_state()
        with patch("manual.squashed_history.MigrationLoader", side_effect=node_not_found(("auth", "0099_old"))):
            with self.assertRaises(NodeNotFoundError):
                history_state()
