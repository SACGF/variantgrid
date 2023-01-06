import json
from abc import abstractmethod, ABC
from dataclasses import dataclass
from decimal import Decimal
from functools import cached_property
from typing import Union, Dict, Any, List, Mapping


# Inclusion of this code snippet will cause "to_json()" to be called on classes by the JSONEncoder, allowing them to become serializable
# START SNIPPET - JSON INCLUSION
def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _default
# END SNIPPET

# JSON Types, just useful for documenting
JsonPrimitiveType = Union[str, int, float, bool, None]
JsonObjType = Dict[JsonPrimitiveType, 'JsonDataType']
JsonListType = List['JSonDataType']
JsonDataType = Union[JsonListType, JsonObjType, JsonPrimitiveType]


def make_json_safe_in_place(obj):
    """
    converts all Decimals in a nested dict/array to floats.
    Specifically because ijson may return some
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, Decimal):
                obj[key] = float(value)
            else:
                make_json_safe_in_place(value)
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            if isinstance(value, Decimal):
                obj[index] = float(value)
            else:
                make_json_safe_in_place(value)
    else:
        pass


def force_json(obj: Any) -> JsonDataType:
    """
    converts an object into a JSONDict - assuming it has the appropriate methods
    does this by converting it to a JSON string and then parsing that string
    not overly efficient but the only way I know how when it comes to nested objects
    that have a JSON representation
    """
    return json.loads(json.dumps(obj))


def strip_json(json_values: JsonDataType) -> JsonDataType:
    """
    Remove null, empty strings and false and empty lists from JSON values
    (contents of arrays won't be affected).
    Can optimise exports of json data when fields are often blank
    """
    if isinstance(json_values, Mapping):
        ret_value = {}
        for key, value in json_values.items():
            if value == '' or value is None or value is False or (isinstance(value, list) and not value):
                pass
            else:
                value = strip_json(value)
                ret_value[key] = value
        return ret_value

    if isinstance(json_values, list):
        ret_value = []
        for value in json_values:
            ret_value.append(strip_json(value))
        return ret_value
    return json_values


class JsonPathPart(ABC):
    @property
    @abstractmethod
    def short(self) -> str:
        pass


@dataclass(frozen=True)
class JsonPathIndex(JsonPathPart):
    index: int

    @property
    def short(self):
        return str(self.index)

    def __str__(self):
        return f'[{self.index}]'


@dataclass(frozen=True)
class JsonPathKey(JsonPathPart):
    key: str

    @property
    def short(self):
        return str(self.key)

    def __str__(self):
        return f"[{json.dumps(self.key)}]"

@dataclass(frozen=True)
class JsonPathId(JsonPathPart):
    key: str
    # not a real path, indicates that we're treating a list as
    # a dictionary with "id" as the key
    @property
    def short(self):
        return f"id({self.key})"

    def __str__(self):
        return f"id({self.key})"


@dataclass(frozen=True, repr=False)
class JsonDiff:
    json_path: List[JsonPathPart]
    a: JsonDataType
    b: JsonDataType

    def __repr__(self) -> str:
        return f"root{self.json_path_str} {self.a} -> {self.b}"

    @cached_property
    def json_path_str(self):
        return "".join(str(p) for p in self.json_path)

    @property
    def json_path_short(self):
        return ".".join(p.short for p in self.json_path)

    def __lt__(self, other):
        # TODO handle proper index differences
        return self.json_path_str < other.json_path_str


class JsonDiffs:

    def __init__(self, json_diffs: List[JsonDiff]):
        self.json_diffs = json_diffs

    def __bool__(self):
        return bool(self.json_diffs)

    def to_json(self, before_label: str = "before", after_label: str = "after") -> JsonObjType:
        diff_dict = {}
        for diff in self.json_diffs:
            diff_dict[diff.json_path_short] = {
                before_label: diff.a,
                after_label: diff.b
            }
        return diff_dict

    @staticmethod
    def differences(obj1: JsonDataType, obj2: JsonDataType) -> 'JsonDiffs':
        diffs: List['JsonDiff'] = []
        JsonDiffs._differences(obj1, obj2, [], diffs)
        diffs.sort()
        return JsonDiffs(diffs)

    @staticmethod
    def _differences(obj1: JsonDataType, obj2: JsonDataType, path: List[JsonPathPart], diffs: List[JsonDiff]) -> None:
        if obj1 == obj2:
            return

        if type(obj1) is type(obj2):
            if isinstance(obj1, dict):
                for key in obj1.keys() | obj2.keys():
                    JsonDiffs._differences(obj1.get(key), obj2.get(key), path + [JsonPathKey(key)], diffs)
            elif isinstance(obj1, list):
                has_id = all(isinstance(x, dict) and x.get('id') is not None for x in obj1 + obj2)
                if has_id:
                    obj1dict = {x.get('id'): x for x in obj1}
                    obj2dict = {x.get('id'): x for x in obj2}
                    for key in obj1dict.keys() | obj2dict.keys():
                        JsonDiffs._differences(obj1dict.get(key), obj2dict.get(key), path + [JsonPathId(key)], diffs)
                else:
                    max_length = max(len(obj1), len(obj2))
                    for index in range(0, max_length):
                        obj1_index = obj1[index] if len(obj1) > index else None
                        obj2_index = obj2[index] if len(obj2) > index else None
                        JsonDiffs._differences(obj1_index, obj2_index, path + [JsonPathIndex(index)], diffs)
            else:
                diffs.append(JsonDiff(json_path=path, a=obj1, b=obj2))
        else:
            diffs.append(JsonDiff(json_path=path, a=obj1, b=obj2))
