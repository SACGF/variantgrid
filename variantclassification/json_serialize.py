import collections
import json


# A little hack that will make call 'to_json' on any non json serializable class
# So you can implement to_json to become json serializable
def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)


_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _default


def strip_json(json_values):
    """
    Remove null, empty strings and false and empty lists from JSON values
    (contents of arrays wont be affected)
    """
    if isinstance(json_values, collections.Mapping):
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
