import logging
from collections import Counter

from django.core.management import BaseCommand

from classification.enums import SpecialEKeys
from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):
        user = admin_bot()

        count = 0
        lab_changes = Counter()
        for classification in Classification.objects.filter(evidence__allele_frequency__isnull=False):
            old_value = classification.get("allele_frequency")
            if old_value is not None:
                # Get out dict - so we can look at and store notes
                value_obj = dict(classification.evidence.get("allele_frequency"))
                existing_note = value_obj.get("note")
                if existing_note:
                    if "Converted from" in existing_note:
                        continue  # Already run

                try:
                    to_value = float(old_value) / 100
                except ValueError:
                    # Someone had entered "0..2"
                    # Someone entered '34%'
                    value = old_value.replace("..", ".").replace("%", "")
                    try:
                        to_value = float(value) / 100
                    except ValueError:
                        logging.error("Couldn't convert classification: %d", classification.pk)
                        raise

                value_obj["value"] = to_value
                notes = []
                existing_note = value_obj.get("note")
                if existing_note:
                    notes.append(existing_note)
                notes.append(f"Converted from '{old_value}'%")
                value_obj["note"] = "\n".join(notes)
                patch = {
                    SpecialEKeys.ALLELE_FREQUENCY: value_obj,
                }
                classification.revalidate(user, migration_patch=patch)

                lab_changes[classification.lab.name] += 1
                count += 1

                if count % 100 == 0:
                    print(f"Processed {count} records")

        print("Classifications changed per lab:")
        print(lab_changes)
