from dataclasses import dataclass
from typing import Any, Optional

from django.core.management import BaseCommand
import re

from classification.models import Classification
from snpdb.models import Lab


@dataclass
class DataFix:
    field: str
    suspect_value: str

    def __str__(self):
        return f"{self.field}: !{self.suspect_value}!"

@dataclass
class DataFixRun:
    type: str
    fixes: list[DataFix]

    def __str__(self):
        return f"{self.type} - " + ", ".join(str(df) for df in self.fixes)


class DataFixer:

    def __init__(self, bad_pattern: re.Pattern, replacement: Any):
        self.bad_pattern = bad_pattern
        self.replacement = replacement
        print(self.replacement)

    def fix_classification_data(self, type: str, data: dict) -> Optional[DataFixRun]:
        modifications = []
        for key, blob in data.items():
            if isinstance(blob, dict):
                if text := blob.get("value"):
                    if isinstance(text, str):
                        if matches := self.bad_pattern.findall(text):
                            modifications.append(DataFix(key, ", ".join(matches)))
                            blob["value"] = self.replacement

        if modifications:
            return DataFixRun(type=type, fixes=modifications)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--lab', type=str, required=True)
        parser.add_argument('--pattern', type=str, required=True)
        parser.add_argument('--pattern_icase', action='store_true')
        parser.add_argument('--apply', type=int, default=0)
        parser.add_argument('--replacement', type=str, default='REDACTED')

    def handle(self, *args, **options):
        lab_id = options["lab"]
        pattern = options["pattern"]
        pattern_icase = options["pattern_icase"]
        apply_remaining = options["apply"]
        replacement = options["replacement"]

        print(f"Lab = {str(Lab.objects.get(id=lab_id))}")
        print(f"Pattern = {pattern}, icase = {pattern_icase}")
        print(f"Apply to = {apply_remaining}")
        print(f"Replacement = \"{replacement}\"")

        bad_pattern = re.compile(pattern, flags=re.IGNORECASE if pattern_icase else 0)
        data_fixer = DataFixer(bad_pattern, replacement)

        for classification in Classification.objects.filter(lab_id=lab_id).iterator():
            changes = []
            if class_change := data_fixer.fix_classification_data("classification", classification.evidence):
                changes.append(class_change)
                if apply_remaining:
                    classification.save(update_fields=["evidence"])

            for cm in classification.classificationmodification_set.order_by('created').all():
                modifying_mod = False
                if published := cm.published_evidence:
                    if published_change := data_fixer.fix_classification_data(f"published {cm.pk}", published):
                        changes.append(published_change)
                        modifying_mod = True
                if delta_changes := data_fixer.fix_classification_data(f"delta {cm.pk}", cm.delta):
                    changes.append(delta_changes)
                    modifying_mod = True

                if modifying_mod and apply_remaining:
                    cm.save(update_fields=["published", "delta"])

            if changes:
                print(f"{classification.friendly_label}: ")
                for change in changes:
                    print("  " + str(change))
                if apply_remaining:
                    print("UPDATED")
                    apply_remaining -= 1