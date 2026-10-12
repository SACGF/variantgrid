from django.core.management import BaseCommand


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        # No longer relevant
        pass
