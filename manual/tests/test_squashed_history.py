from unittest.mock import patch

from django.core.management import CommandError, call_command
from django.db import connection
from django.db.migrations.recorder import MigrationRecorder
from django.test import TestCase

from manual.squashed_history import (
    CHECKPOINT_MIGRATION,
    HistoryState,
    history_state,
    refuse_celery_worker_on_old_history,
    reset_squashed_history,
)

OLD_MIGRATION = ("snpdb", "0297_retire_somalier_download_manual_task")


class SquashedHistoryTest(TestCase):

    def setUp(self):
        self.recorder = MigrationRecorder(connection)

    def test_fresh_database_is_current(self):
        self.assertEqual(history_state(), HistoryState.CURRENT)

    def test_old_history_without_checkpoint_refuses_reset(self):
        self.recorder.record_applied(*OLD_MIGRATION)
        self.assertEqual(history_state(), HistoryState.NEEDS_CHECKPOINT)
        with self.assertRaises(CommandError):
            call_command("reset_squashed_migration_history")

    def test_reset_swaps_old_rows_for_initials(self):
        self.recorder.record_applied(*OLD_MIGRATION)
        self.recorder.record_applied(*CHECKPOINT_MIGRATION)
        self.recorder.migration_qs.filter(app="snpdb", name="0001_initial").delete()
        self.assertEqual(history_state(), HistoryState.NEEDS_RESET)

        reset_squashed_history()

        self.assertEqual(history_state(), HistoryState.CURRENT)
        applied = self.recorder.applied_migrations()
        self.assertNotIn(OLD_MIGRATION, applied)
        self.assertIn(("snpdb", "0001_initial"), applied)
        self.assertIn(CHECKPOINT_MIGRATION, applied)

    def test_celery_worker_exits_on_old_history(self):
        """ celery swallows an Exception from a worker_init handler, so only SystemExit stops the worker """
        with patch("manual.squashed_history.refuse_to_serve_old_history", side_effect=RuntimeError("old history")):
            with self.assertRaises(SystemExit):
                refuse_celery_worker_on_old_history()
