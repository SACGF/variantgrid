from datetime import timedelta

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from snpdb.models.models_api_token import ApiToken


class Command(BaseCommand):
    category = "ops"
    help = "Create an API token for a user and print its key - the way to issue one where users can't make their own"

    def add_arguments(self, parser):
        parser.add_argument('--username', required=True)
        parser.add_argument('--name', required=True, help="What the token is for, shown on the API token pages")
        parser.add_argument('--days', type=int, help="Days until it expires (default: never)")

    def handle(self, *args, **options):
        user = User.objects.get(username=options["username"])
        expiry = timedelta(days=days) if (days := options["days"]) else None
        _token, key = ApiToken.create_for_user(user, options["name"], expiry, created_by=None)
        print(key)
