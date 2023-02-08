from typing import List, Optional, Dict, Any
from django.http import HttpRequest
from classification.models import EvidenceKey, EvidenceKeyMap
from classification.models.evidence_mixin import VCBlobDict
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter2
import json

"""
For generating a report about the usage of evidence keys.
There is code that keeps track of the first record id for each value, but for the sake of simplicity
the report just prints how many time each value was encountered.
"""

MAX_EXAMPLES = 15
MAX_EXAMPLE_LEN = 50


class ValueExample:

    def __init__(self, value: Any, example_source: Optional[Any]):
        self.value = value
        self.count = 1
        self.example_source = example_source

    def to_json(self) -> Dict:
        data = {"value": self.value, "count": self.count}
        return data


class ValueCounter:

    def __init__(self):
        self.used = 0
        self.unused = 0
        self.others = 0
        self.examples: Optional[Dict[str, ValueExample]] = None

    def _count_example(self, value: Any, example_source: Optional[Any] = None, from_part: bool = False):
        key = str(value)
        existing = self.examples.get(key)
        if existing:
            existing.count += 1
        else:
            existing = ValueExample(value, example_source)
            self.examples[key] = existing
        if from_part:
            existing.from_part = (existing.from_part if existing.from_part else 0) + 1

    def count(self, value: Any, source: Optional[Any] = None):
        if value is None:
            self.unused += 1
            return

        self.used += 1
        if self.examples is None:
            self.examples = {}

        if isinstance(value, list):
            """
            if len(value) > 1:
                for sub_value in value:
                    self._count_example(value = sub_value, example_source=source, from_part = True)
                self._count_example(value = value, example_source=source)
            elif len(value) == 1:
                self._count_example(value = value, example_source=source)
            """
            value = ", ".join([str(item) for item in value])
            self._count_example(value=value, example_source=source)

        elif isinstance(value, bool):
            self._count_example(value=value, example_source=source)

        else:
            if isinstance(value, str) and len(value) > MAX_EXAMPLE_LEN:
                value = value[0:MAX_EXAMPLE_LEN] + '...(' + str(len(value)) + ')chars'

            has_match = value in self.examples
            if has_match or len(self.examples) < MAX_EXAMPLES:
                self._count_example(value=value, example_source=source)
            else:
                self.others += 1

    def to_json(self):
        data: Dict = {}
        if self.used:
            percent = (self.used / (self.used + self.unused)) * 100
            data["used_percentage"] = f"{percent:.2f}%"
        else:
            return "0%"
        if self.examples:
            example_list = [example for example in self.examples.values()]
            example_list.sort(key=lambda x: x.count, reverse=True)
            """
            data["examples"] = [example.to_json() for example in example_list]
            """
            example_dict = {}
            for example in example_list:
                example_dict[str(example.value)] = example.count
            if self.others:
                example_dict["<others>"] = self.others
            data["value_counts"] = example_dict

        return data


class KeyCount:

    def __init__(self, e_key: EvidenceKey):
        self.e_key = e_key
        self.values = ValueCounter()
        self.notes = ValueCounter()
        self.explain = ValueCounter()
        self.dbrefs: Dict[str, int] = {}

    def count(self, blob: Optional[VCBlobDict], source: Any):
        if not blob:
            self.values.count(value=None)
            self.notes.count(value=None)
            self.explain.count(value=None)
        else:
            self.values.count(value=blob.get('value'), source=source)
            self.notes.count(value=blob.get('note'), source=source)
            self.explain.count(value=blob.get('explain'), source=source)
            if refs := blob.get('db_refs'):
                for ref in refs:
                    db = ref.get('db')
                    self.dbrefs[db] = self.dbrefs.get(db, 0) + 1

    def to_json(self):
        if not self.values.used and not self.notes.used and not self.explain.used:
            return False
        data = {
            "values": self.values.to_json()
        }
        if self.notes.used:
            data["notes"] = self.notes.to_json()
        if self.explain.used:
            data["explains"] = self.explain.to_json()
        if self.dbrefs:
            data['reference_counts'] = self.dbrefs
        return data

@register_classification_exporter("keys")
class ClassificationExportFormatterKeys(ClassificationExportFormatter2):
    """
    Exports data in the format that Agilent's Alissa can import it
    """

    def __init__(self, classification_filter: ClassificationFilter):
        key_counters: Dict[str, KeyCount] = {}
        for e_key in EvidenceKeyMap.instance().all_keys:
            key_counters[e_key.key] = KeyCount(e_key=e_key)
        self.key_counters = key_counters
        self.rows_inspected = 0
        self.rows_with_errors = 0
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterKeys':
        classification_filter = ClassificationFilter.from_request(request)
        return ClassificationExportFormatterKeys(
            classification_filter=classification_filter
        )

    def header(self) -> List[str]:
        return None

    def row(self, data: AlleleData) -> List[str]:
        self.rows_with_errors += len([issue for issue in data.all_cms if issue.validation_include])
        for ci in data.all_cms:
            if ci.withdrawn:
                continue # don't count withdrawn
            if ci.has_issue:
                # still count towards total
                self.rows_with_errors += 1

            self.rows_inspected += 1
            classification = ci.classification
            data = classification.evidence
            for key, counter in self.key_counters.items():
                blob = data.get(key)
                counter.count(blob, source=classification.id)

    def footer(self) -> List[str]:
        self.row_count = 1
        data = {}
        for key, value in self.key_counters.items():
            data[key] = value.to_json()

        keys = {
            "classifications_inspected": self.rows_inspected,
            "classifications_included_but_excluded_from_download": self.rows_with_errors,
            "keys": data
        }

        return [json.dumps(keys, indent=4)]

    def content_type(self) -> str:
        return "application/json"

    def extension(self) -> str:
        return "json"
