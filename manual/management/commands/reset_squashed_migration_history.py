from django.core.management.base import BaseCommand, CommandError

from manual.squashed_history import (
    HistoryState,
    history_state,
    history_state_message,
    reset_squashed_history,
)


class Command(BaseCommand):
    """ Once per deployment, after branch migration_checkpoint: swaps the pre-squash migration rows for the new
        initials (#736). The upgrader runs it before migrate; a no-op on a current database """
    category = "ops"
    requires_system_checks = []  # check_squashed_history errors on exactly the database this command fixes

    def handle(self, *args, **options):
        state = history_state()
        if state == HistoryState.CURRENT:
            self.stdout.write(history_state_message(state))
        elif state == HistoryState.NEEDS_CHECKPOINT:
            raise CommandError(history_state_message(state))
        else:
            deleted = reset_squashed_history()
            self.stdout.write(f"Deleted {deleted} pre-squash migration rows and recorded the new initial migrations")
