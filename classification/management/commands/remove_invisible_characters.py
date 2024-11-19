import json

from django.core.management.base import BaseCommand
import re
from classification.models import Classification
from library.guardian_utils import admin_bot


def has_invisible_characters(text):
    pattern = "\u00a0"
    return re.search(pattern, text)


def ensure_string(data):
    if isinstance(data, (dict, list)):
        return json.dumps(data)
    elif isinstance(data, str):
        return data
    else:
        return str(data)


class Command(BaseCommand):

    def handle(self, *args, **options):
        classifications = Classification.objects.all()
        user = admin_bot()

        for classification in classifications:
            evidence = classification.evidence

            for key, value in evidence.items():
                for k, v in value.items():
                    if match := has_invisible_characters(ensure_string(v)):
                        print(f"match found in {classification.id} record, in: {key}")
                        classification.revalidate(user=user)

        self.stdout.write(
            self.style.SUCCESS('Invisible characters removed from evidence field for all classifications.'))
