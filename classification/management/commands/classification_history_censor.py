from typing import Any

from django.core.management import BaseCommand
import re

from classification.models import Classification
from snpdb.models import Lab


REPLACEMENT_REDACT = object()

class DataFixer:

    def __init__(self, bad_pattern: re.Pattern, replacement: Any):
        self.bad_pattern = bad_pattern
        self.replacement = replacement

    def fix_classificaiton_data(self, data: dict) -> bool:
        has_modified = False
        for key, blob in data.items():
            if isinstance(blob, dict):
                if text := blob.get("value"):
                    if isinstance(text, str):
                        if self.bad_pattern.search(text):
                            has_modified = True
                            if self.replacement == REPLACEMENT_REDACT:
                                blob["value"] = "REDACTED"
                            else:
                                raise ValueError(f"Replacement method {self.replacement} is not supported")
        return has_modified


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
            classification_fix = False
            history_count = 0
            if data_fixer.fix_classificaiton_data(classification.evidence):
                classification_fix = True

            for cm in classification.classificationmodification_set.all():
                modification_fix = False
                if published := cm.published_evidence:
                    if data_fixer.fix_classificaiton_data(published):
                        modification_fix = True
                    if data_fixer.fix_classificaiton_data(cm.delta):
                        modification_fix = True
                if modification_fix:
                    history_count += 1

            if classification_fix or history_count:
                print(f"{classification.friendly_label} Classification fixed = {classification_fix}, history fixes = {history_count}")
