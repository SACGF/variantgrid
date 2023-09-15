from django.core.management.base import BaseCommand

from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):
        user = admin_bot()
        classifications = Classification.objects.filter(
            condition_resolution__isnull=False,
            evidence__condition__validation__icontains="error"
        )
        print(f"Found {classifications.count()} in database")
        for classification in classifications:
            print(f"Revalidating {classification}")
            classification.revalidate(user=user)
        self.stdout.write(self.style.SUCCESS("Revalidation completed to override blank condition."))
