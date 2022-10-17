from decimal import Decimal
from typing import Union, Dict, Any, List

# this needs to be moved out of uicore and into json_utils
from uicore.json.json_types import JsonObjType

# isn't restrictive enough, but helps documention
# DEPRECATED use JsonObjType or similar
JSON = Union[Dict[str, Any], List[Any], int, float, str, bool]


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
