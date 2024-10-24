from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from rest_framework.authtoken.models import Token


class Command(BaseCommand):
    help = 'Copy data from a one VCF to another (with the same samples)'

    def add_arguments(self, parser):
        parser.add_argument('--username', required=True)

    def handle(self, *args, **options):
        username = options["username"]
        user = User.objects.get(username=username)
        token = Token.objects.create(user=user)
        print(token.key)