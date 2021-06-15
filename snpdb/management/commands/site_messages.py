import datetime
from argparse import ArgumentError

from django.core.management.base import BaseCommand
from django.utils import timezone

from library.enums.log_level import LogLevel
from library.utils import invert_dict
from snpdb.models import SiteMessage


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--shutdown', type=int, required=False)
        parser.add_argument('--message', required=False)
        parser.add_argument('--severity', help='DEBUG/INFO/WARNING/ERROR', default='INFO')
        parser.add_argument('--clear-old', action='store_true', help='Clear messages with date_time in the past')
        parser.add_argument('--clear-all', action='store_true')

    def handle(self, *args, **options):
        clear_all = options.get("clear_all")
        clear_old = options.get("clear_old")
        shutdown = options.get("shutdown")
        message = options.get("message")
        severity = options["severity"]
        log_level_choices = invert_dict(dict(LogLevel.CHOICES))
        log_level = log_level_choices.get(severity)
        if log_level is None:
            raise ArgumentError(f"Severity must be one of '{','.join(log_level_choices.keys())}'")

        if clear_all:
            print("Deleting all site messages")
            SiteMessage.objects.all().delete()
        elif clear_old:
            for sm in SiteMessage.objects.filter(date_time__lt=timezone.now()):
                print(f"Deleting old site message: '{sm}'")
                sm.delete()

        if shutdown:
            date_time = timezone.now() + datetime.timedelta(minutes=shutdown)
            sm = SiteMessage.objects.create(message="The system will be shut down for maintenance",
                                            log_level=LogLevel.ERROR,
                                            date_time=date_time)
            print(f"Added message: '{sm}'")

        if message:
            sm = SiteMessage.objects.create(message=message, log_level=log_level)
            print(f"Added message: '{sm}'")
