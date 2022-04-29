import json
# Inclusion of this code snippet will cause "to_json()" to be called on classes by the JSONEncoder, allowing them to become serializable
from typing import Any

from uicore.json.json_types import JsonDataType


def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _default


def force_json(obj: Any) -> JsonDataType:
    # converts an object into a JSONDict - assuming it has the appropriate methods
    # does this by converting it toa JSON string and then parsing that string
    # not overy efficient but the only way I know how when it comes to nested objects
    # that have a JSON representation
    return json.loads(json.dumps(obj))