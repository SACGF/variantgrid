import socket
from dataclasses import dataclass
import datetime
from typing import Any, Optional

from django.core.management import BaseCommand
import re

from classification.models import Classification
from snpdb.models import Lab, Organization

VERSION = 2

@dataclass
class DataFix:
    field: str
    suspect_values: list[str]

    def __str__(self):
        return f"{self.field}: {', '.join(self.suspect_values)}"

@dataclass
class DataFixRun:
    type: str
    fixes: list[DataFix]

    def __str__(self):
        return f"{self.type} - " + ", ".join(str(df) for df in self.fixes)


class DataFixer:

    def __init__(self, bad_pattern: re.Pattern, replacement: Any, keys: Optional[set[str]] = None):
        self.bad_pattern = bad_pattern
        self.replacement = replacement
        self.keys = keys

    def fix_classification_data(self, type: str, data: dict) -> Optional[DataFixRun]:
        modifications = []
        for key, blob in data.items():
            if self.keys and key not in self.keys:
                continue
            if isinstance(blob, dict):
                if text := blob.get("value"):
                    if isinstance(text, str):
                        suspect_values = []
                        for match in self.bad_pattern.finditer(text):
                            start = match.start()
                            while start > 0 and text[start-1] not in {" ", "\t", "\n", ">"}:
                                start -= 1
                            end = match.end()
                            while end < len(text) - 1 and text[end] not in {" ", "\t", "\n", ">"}:
                                end += 1
                            suspect_value = text[start:end]
                            suspect_values.append(suspect_value)

                        if suspect_values:
                            modifications.append(DataFix(field=key, suspect_values=suspect_values))
                            blob["value"] = self.replacement


        if modifications:
            return DataFixRun(type=type, fixes=modifications)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--lab', type=str, required=False)
        parser.add_argument('--org', type=str, required=False)
        parser.add_argument('--pattern', type=str, required=True)
        parser.add_argument('--pattern_icase', action='store_true')
        parser.add_argument('--apply', type=int, default=0)
        parser.add_argument('--keys', type=str, default=None)
        parser.add_argument('--replacement', type=str, required=True)

    def handle(self, *args, **options):
        lab_id:str = options["lab"]
        org_id: str = options["org"]
        pattern = options["pattern"]
        pattern_icase = options["pattern_icase"]
        apply_remaining = options["apply"]
        replacement = options["replacement"]

        print(f"History Redactor v {VERSION}")
        print(f"Server: {socket.gethostname()}")
        print(f"Server Time: {datetime.datetime.now()}")

        lab: Optional[Lab] = None
        org: Optional[Organization] = None

        if lab_id:
            if lab_id.isnumeric():
                lab = Lab.objects.get(id=lab_id)
            else:
                lab = Lab.objects.get(group_name=lab_id)
            print(f"Lab = {str(lab)}")

        if org_id:
            if org_id.isnumeric():
                org = Organization.objects.get(id=org_id)
            else:
                org = Organization.objects.get(group_name=org_id)
            print(f"Org = {str(org)}")


        print(f"Pattern = {pattern}, icase = {pattern_icase}")
        if not apply_remaining:
            print("Dry Run - NO UPDATES")
        else:
            print(f"Apply to = {apply_remaining}")
        print(f"Replacement = \"{replacement}\"")

        print("")
        keys = None
        if keys_str := options["keys"]:
            keys = set([key.strip() for key in keys_str.split(",")])

        bad_pattern = re.compile(pattern, flags=re.IGNORECASE if pattern_icase else 0)
        data_fixer = DataFixer(bad_pattern, replacement, keys=keys)

        qs = Classification.objects.all()
        if lab:
            qs = qs.filter(lab=lab)
        if org:
            qs = qs.filter(lab__organization=org)

        for classification in qs.iterator():
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
                    cm.save(update_fields=["published_evidence", "delta"])

            if changes:
                print(f"{classification.friendly_label}: ")
                for change in changes:
                    print("  " + str(change))
                if apply_remaining:
                    print("UPDATED")
                    apply_remaining -= 1

        print("FINISHED")