from collections import Counter
from typing import Mapping

from django.core.management import BaseCommand

from classification.enums import SubmissionSource, SpecialEKeys
from classification.models import Classification
from library.guardian_utils import admin_bot


class Command(BaseCommand):

    def handle(self, *args, **options):
        user = admin_bot()

        lab_changes = Counter()
        for classification in Classification.objects.filter(evidence__allele_frequency__isnull=False):
            old_value = classification.get("allele_frequency")
            if old_value is not None:
                # Get out dict - so we can look at and store notes
                value_obj = classification.evidence.get("allele_frequency")
                if not isinstance(value_obj, Mapping):
                    value_obj = {}
                # Someone had entered "0..2"
                old_value = old_value.replace("..", ".")
                to_value = str(float(old_value) / 100)
                value_obj["value"] = to_value
                notes = []
                existing_note = value_obj.get("note")
                if existing_note:
                    if "Converted from" in existing_note:
                        raise ValueError(f"Classification {classification} contains allele_frequency note with "
                                         f"'converted from' - Legacy upgrade appears to have been run twice!!")
                    notes.append(existing_note)
                notes.append(f"Converted from '{old_value}'%")
                value_obj["note"] = ". ".join(notes)
                patch = {
                    SpecialEKeys.ALLELE_FREQUENCY: value_obj,
                }
                classification.patch_value(patch=patch,
                                           source=SubmissionSource.VARIANT_GRID,
                                           save=True,
                                           user=user)
                classification.revalidate(user)
                lab_changes[classification.lab.name] += 1

        print("Classifications changed per lab:")
        print(lab_changes)
