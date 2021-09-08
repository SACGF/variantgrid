from dateutil.relativedelta import relativedelta
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from django.utils import timezone


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--weeks', type=int, default=26)
        parser.add_argument('--separator', default=";")

    def handle(self, *args, **options):
        weeks = options.get("weeks")
        separator = options.get("separator")

        last_login_date = timezone.now() + relativedelta(weeks=-weeks)
        print(f"Users logged in during last {weeks} weeks (after {last_login_date})")

        users_without_emails = []
        emails = set()

        for user in User.objects.filter(last_login__gte=last_login_date):
            if user.email:
                emails.add(user.email)
            else:
                users_without_emails.append(user.username)

        if users_without_emails:
            print("Users without emails:")
            print(", ".join(users_without_emails))

        if emails:
            print("Emails:")
            print(separator.join(emails))
