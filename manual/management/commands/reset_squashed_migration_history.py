from django.core.management.base import BaseCommand, CommandError

from manual.squashed_history import (
    CheckoutBeforeSquashError,
    HistoryState,
    history_state,
    history_state_message,
    reset_squashed_history,
)


class Command(BaseCommand):
    """ Once per deployment, after branch migration_checkpoint: swaps the pre-squash migration rows for the new
        initials (#736). Run by hand, before scripts/upgrade.sh; a no-op on a current database """
    category = "ops"
    requires_system_checks = []  # check_squashed_history errors on exactly the database this command fixes

    def handle(self, *args, **options):
        try:
            state = history_state()
        except CheckoutBeforeSquashError as e:
            raise CommandError(str(e)) from e
        if state == HistoryState.CURRENT:
            self.stdout.write(history_state_message(state))
        elif state == HistoryState.NEEDS_CHECKPOINT:
            raise CommandError(history_state_message(state))
        else:
            deleted = reset_squashed_history()
            self.stdout.write(f"Deleted {deleted} pre-squash migration rows and recorded the new initial migrations")
