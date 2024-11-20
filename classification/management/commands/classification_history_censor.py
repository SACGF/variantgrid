from dataclasses import dataclass
from typing import Any, Optional

from django.core.management import BaseCommand
import re

from classification.models import Classification
from snpdb.models import Lab


REPLACEMENT_REDACT = object()


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

    def fix_classification_data(self, type: str, data: dict) -> Optional[DataFixRun]:
        modifications = []
        for key, blob in data.items():
            if isinstance(blob, dict):
                if text := blob.get("value"):
                    if isinstance(text, str):
                        if self.bad_pattern.search(text):
                            modifications.append(DataFix(key, ", ".join(self.bad_pattern.findall(text))))
                            if self.replacement == REPLACEMENT_REDACT:
                                blob["value"] = "REDACTED"
                            else:
                                raise ValueError(f"Replacement method {self.replacement} is not supported")

        if modifications:
            return DataFixRun(type=type, fixes=modifications)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--lab', type=str, required=True)
        parser.add_argument('--pattern', type=str, required=True)

    def handle(self, *args, **options):
        lab_id = options["lab"]
        pattern = options["pattern"]
        bad_pattern = re.compile(pattern)

        data_fixer = DataFixer(bad_pattern, REPLACEMENT_REDACT)

        for classification in Classification.objects.filter(lab_id=lab_id).iterator():
            changes = []
            if class_change := data_fixer.fix_classification_data("classification", classification.evidence):
                changes.append(class_change)

            for cm in classification.classificationmodification_set.order_by('created').all():
                if published := cm.published_evidence:
                    if published_change := data_fixer.fix_classification_data(f"published {cm.pk}", published):
                        changes.append(published_change)
                    if delta_changes := data_fixer.fix_classification_data(f"delta {cm.pk}", cm.delta):
                        changes.append(delta_changes)

            if changes:
                print(f"{classification.friendly_label}: ")
                for change in changes:
                    print("  " + str(change))