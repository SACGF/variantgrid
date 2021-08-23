import json

# Inclusion of this code snippet will cause "to_json()" to be called on classes by the JSONEncoder, allowing them to become serializable


def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _default
