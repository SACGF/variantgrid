import json
from typing import Dict, Any, Optional

from django.db.models import QuerySet
from django.http import StreamingHttpResponse

from classification.models import ClassificationModification, EvidenceKey, EvidenceKeyMap
from classification.models.evidence_mixin import VCBlobDict
from classification.views.classification_export_utils import BaseExportFormatter
from library.utils import DebugTimer

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
            value = ", ".join(value)
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


class ExportFormatterKeys(BaseExportFormatter):
    """
    Formats as a general CSV
    """

    def benchmark(self, row_limit: int = 100) -> DebugTimer:
        raise NotImplementedError("Benchmark not supported by ExportFormatterKeys")

    def __init__(self, qs: QuerySet, *args, **kwargs):
        self.qs = qs.order_by("-modified")
        self.key_counters: Dict[str, KeyCount] = {}

        super().__init__(*args, **kwargs)

        for e_key in EvidenceKeyMap.instance().all_keys:
            self.key_counters[e_key.key] = KeyCount(e_key=e_key)

    def process(self):
        for vcm in self.qs.iterator(chunk_size=1000):
            self.count_classification(vcm)

    def export(self, as_attachment: bool = True) -> StreamingHttpResponse:
        def stream_response():
            yield '{"keys":'
            self.process()

            data = {}
            for key, value in self.key_counters.items():
                data[key] = value.to_json()

            yield json.dumps(data, indent=4)
            yield '}'
            return None

        response = StreamingHttpResponse(stream_response(), content_type=self.content_type())
        if as_attachment:
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
        return response

    def count_classification(self, classification: ClassificationModification):
        data = classification.evidence
        for key, counter in self.key_counters.items():
            blob = data.get(key)
            counter.count(blob, source=classification.id)

    def content_type(self) -> str:
        return 'application/json'

    def filename(self) -> str:
        return "evidence_key_report.json"
