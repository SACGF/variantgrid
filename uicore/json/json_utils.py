import collections


def strip_json(json_values):
    """
    Remove null, empty strings and false and empty lists from JSON values
    (contents of arrays wont be affected).
    Can optimise exports of json data when fields are often blank
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