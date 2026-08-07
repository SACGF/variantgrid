from django.core.management.base import BaseCommand

from snpdb.models import JobsControl


class Command(BaseCommand):
    help = ("Pause / resume / show the analysis + annotation job dispatchers (operational safety "
            "brake). Pausing stops new work being leased or launched; in-flight work is untouched.")

    def add_arguments(self, parser):
        parser.add_argument("action", choices=["pause", "resume", "status"])
        parser.add_argument("--reason", default="", help="Note shown in status (for 'pause')")

    def handle(self, *args, **options):
        action = options["action"]
        by = f"manage.py jobs_control {action}"
        if action == "pause":
            obj = JobsControl.pause(reason=options["reason"] or "Paused via management command", by=by)
            self.stdout.write(self.style.WARNING(f"PAUSED: {obj.reason}"))
        elif action == "resume":
            JobsControl.resume(by=by)
            self.stdout.write(self.style.SUCCESS("RESUMED - dispatchers will pick work back up"))
        else:  # status
            obj = JobsControl.get_solo()
            if obj.paused:
                self.stdout.write(self.style.WARNING(
                    f"PAUSED at {obj.paused_at} by {obj.paused_by}: {obj.reason}"))
            else:
                self.stdout.write(self.style.SUCCESS("RUNNING"))
