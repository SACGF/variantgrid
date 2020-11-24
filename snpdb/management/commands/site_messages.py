import datetime
from django.core.management.base import BaseCommand
from django.utils import timezone

from snpdb.models import SiteMessage


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--shutdown', type=int, required=False)
        parser.add_argument('--message', required=False)
        parser.add_argument('--clear-old', action='store_true', help='Clear messages with date_time in the past')
        parser.add_argument('--clear-all', action='store_true')

    def handle(self, *args, **options):
        clear_all = options.get("clear_all")
        clear_old = options.get("clear_old")
        shutdown = options.get("shutdown")
        message = options.get("message")

        if clear_all:
            print("Deleting all site messages")
            SiteMessage.objects.all().delete()
        elif clear_old:
            for sm in SiteMessage.objects.filter(date_time__lt=timezone.now()):
                print(f"Deleting old site message: '{sm}'")
                sm.delete()

        if shutdown:
            date_time = timezone.now() + datetime.timedelta(minutes=shutdown)
            sm = SiteMessage.objects.create(message="The system will soon be shut down soon (approx: %(time_away)s)",
                                            date_time=date_time)
            print(f"Added message: '{sm}'")

        if message:
            sm = SiteMessage.objects.create(message=message)
            print(f"Added message: '{sm}'")
