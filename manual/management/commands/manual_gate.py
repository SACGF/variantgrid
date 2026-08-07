from django.core.management import CommandParser
from django.core.management.base import BaseCommand, CommandError

from library.git import Git
from manual.gates import GATES, gate_status
from manual.models import ManualGateSatisfied


class Command(BaseCommand):
    """ List manual-migration gates, or confirm a ManualGate as satisfied.

        Gates are prerequisites that manual `manage` tasks depend on (see manual.gates). Auto
        gates are checked in code; manual gates (e.g. an external data upgrade) must be confirmed
        here once they are done, which unblocks the tasks that require them. """

    def add_arguments(self, parser: CommandParser):
        parser.add_argument('--satisfy', metavar="GATE_NAME",
                            help="Record a (manual) gate as satisfied")
        parser.add_argument('--note', help="Optional note to record against --satisfy")

    def handle(self, *args, **options):
        if name := options.get("satisfy"):
            gate = GATES.get(name)
            if gate is None:
                raise CommandError(f"Unknown gate '{name}'. Known gates: {', '.join(sorted(GATES))}")
            if gate.kind != "manual":
                raise CommandError(f"Gate '{name}' is a {gate.kind} gate - it is checked in code, "
                                   f"not confirmed by an operator.")
            ManualGateSatisfied.objects.create(name=name, note=options.get("note"),
                                               source_version=Git().hash)
            self.stdout.write(f"Gate '{name}' marked satisfied.")
            return

        for name in sorted(GATES):
            status = gate_status(name)
            mark = "OK " if status.satisfied else "-- "
            self.stdout.write(f"{mark}[{status.kind}] {name}: {status.description}")
