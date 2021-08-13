from typing import Union, Dict, Any, List

# JSON Types, more useful just as documenting than for static type checking

JsonPrimitiveType = Union[str, int, float, bool, None]
JsonObjType = Dict[JsonPrimitiveType, 'JsonDataType']
JsonListType = List['JSonDataType']
JsonDataType = Union[JsonListType, JsonObjType, JsonPrimitiveType]