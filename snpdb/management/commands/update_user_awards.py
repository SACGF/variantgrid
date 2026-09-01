from django.core.management.base import BaseCommand

from snpdb.models.models_enums import AwardPeriod
from snpdb.user_award_updates import update_user_awards
from snpdb.user_awards import get_award_definitions


class Command(BaseCommand):
    help = "Recompute user titles and badges (#1819) - the beat task does this hourly/nightly"

    def add_arguments(self, parser):
        parser.add_argument("--periods", nargs="+", choices=[p.value for p in AwardPeriod],
                            default=[p.value for p in AwardPeriod],
                            help="Title periods to recompute (A=all time, M=month, D=day). Default: all")
        parser.add_argument("--no-badges", action="store_true", help="Skip badge recompute")

    def handle(self, *args, **options):
        periods = [AwardPeriod(p) for p in options["periods"]]
        definitions = get_award_definitions()
        if update_user_awards(periods, badges=not options["no_badges"]):
            self.stdout.write(f"Updated {len(definitions)} award definition(s) for periods {[p.label for p in periods]}")
        else:
            self.stdout.write("USER_AWARDS_ENABLED is off - nothing to do")
