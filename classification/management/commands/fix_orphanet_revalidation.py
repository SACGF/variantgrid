from django.core.management.base import BaseCommand
from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):
        user = admin_bot()
        search_term_orpha = "ORPHANET"
        classifications = Classification.objects.filter(evidence__search_terms__icontains=search_term_orpha)
        print(f"Found {classifications.count()} classifications with '{search_term_orpha}' in literature")
        for classification in classifications:
            classification.revalidate(user=user)
        self.stdout.write(self.style.SUCCESS("Revalidation complete for Orphanet classifications."))
