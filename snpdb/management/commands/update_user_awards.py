from django.core.management.base import BaseCommand

from snpdb.user_award_updates import update_user_awards
from snpdb.user_awards import get_award_definitions


class Command(BaseCommand):
    category = "ops"
    help = "Recompute user badges (#1819) - the beat task does this nightly"

    def handle(self, *args, **options):
        if update_user_awards():
            self.stdout.write(f"Updated {len(get_award_definitions())} award definition(s)")
        else:
            self.stdout.write("USER_AWARDS_ENABLED is off - nothing to do")
