from django.core.management import BaseCommand

from classification.enums import SubmissionSource, SpecialEKeys
from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):
        user = admin_bot()

        modified_classifications = []
        for classification in Classification.objects.filter(evidence__variant_coordinate__value__icontains=".."):
            modified_classifications.append(str(classification.pk))
            patch = {
                SpecialEKeys.VARIANT_COORDINATE: classification.variant.full_string
            }
            classification.patch_value(patch=patch,
                                       source=SubmissionSource.VARIANT_GRID,
                                       save=True,
                                       user=user)
            classification.revalidate(user)

        print(f"Modified {len(modified_classifications)} classifications: {','.join(modified_classifications)}")
