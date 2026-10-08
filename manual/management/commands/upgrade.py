import argparse
import sys

from django.core.management import CommandParser
from django.core.management.base import BaseCommand

from manual.upgrader import Upgrader


class Command(BaseCommand):
    """ The deploy upgrader - run it through ./scripts/upgrade.sh, which installs requirements first """
    category = "ops"
    help = "Deploy upgrader: standard steps (git pull, migrate, collectstatic, deployment_check, deployed) " \
           "and the outstanding manual steps migrations registered, all in this one process. " \
           "With no option, an interactive menu. Set VG_INSTALL_REQUIREMENTS=0 to skip installing requirements."

    def add_arguments(self, parser: CommandParser):
        mode = parser.add_mutually_exclusive_group()
        mode.add_argument("--quick", action="store_true",
                          help="Run the standard steps, stopping at the first failure, and quit if nothing else is "
                               "outstanding (otherwise the menu)")
        mode.add_argument("--auto-manage", action="store_true",
                          help="Run every unblocked 'manage' step, re-evaluating between passes so ordering gates "
                               "release as prerequisites finish, carrying on past failures. Gated, obsolete and "
                               "non-manage steps are reported, not run")
        mode.add_argument("--steps", metavar="KEYS",
                          help="Run these menu keys in order, carrying on past failures (skipping any step gated "
                               "after one that failed), e.g. 'm,c,2-5'")
        mode.add_argument("--resume", help=argparse.SUPPRESS)  # the rest of a selection, after a pull restarted us

    def handle(self, *args, **options):
        upgrader = Upgrader()
        try:
            if options["quick"]:
                exit_code = upgrader.quick()
            elif options["auto_manage"]:
                exit_code = 1 if upgrader.auto_manage() else 0
            elif steps := options["steps"]:
                exit_code = upgrader.run_steps(steps)
            elif resume := options["resume"]:
                exit_code = upgrader.resume(resume)
            else:
                exit_code = upgrader.prompt()
        except KeyboardInterrupt:
            exit_code = 130
        sys.exit(exit_code)
