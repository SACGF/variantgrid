from django.core.management import BaseCommand


class Command(BaseCommand):

    category = "maintenance"

    def handle(self, *args, **options):
        # No longer relevant
        pass
