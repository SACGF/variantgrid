import copy
import json
from dataclasses import dataclass, field
from typing import List, Iterator, Union, Optional

from lazy import lazy

from uicore.json.json_types import JsonDataType


# ValidatedJson (with JSonMessages) is used for serializing to JSon where there's also the need to providing infos or warnings
# e.g. for converting a classification to the ClinVar format, but providing errors around illegal values.
# The "pure" JSON can be extracted, or the ValidatedJson can be serialized (and later rendered)


@dataclass(frozen=True)
class JsonMessage:

    severity: str
    text: str

    @property
    def is_error(self) -> bool:
        # TODO get rid of all these hardcoded strings
        return self.severity == "error"

    @property
    def bs(self) -> str:
        if self.severity == "error":
            return "danger"
        if self.severity == "warning":
            return "warning"
        if self.severity == "info":
            return "info"
        return "info"

    def to_json(self):
        return {"severity": self.severity, "text": self.text}

    @staticmethod
    def deserialize(json_data: JsonDataType):
        return JsonMessage(severity=json_data.get("severity"), text=json_data.get("text"))


@dataclass(frozen=True)
class JsonMessages:

    messages: List[JsonMessage] = field(default_factory=list)

    def errors(self) -> List[JsonMessage]:
        return JsonMessages(messages=[mess for mess in self.messages if mess.severity == "error"])

    def warnings(self) -> List[JsonMessage]:
        return JsonMessages(messages=[mess for mess in self.messages if mess.severity == "warning"])

    @staticmethod
    def error(message: str):
        return JsonMessages([JsonMessage(severity="error", text=message)])

    @staticmethod
    def warning(message: str):
        return JsonMessages([JsonMessage(severity="warning", text=message)])

    @staticmethod
    def info(message: str):
        return JsonMessages([JsonMessage(severity="info", text=message)])

    @staticmethod
    def severity(severity: str, message: str):
        return JsonMessages([JsonMessage(severity=severity, text=message)])

    def __add__(self, other) -> 'JsonMessages':
        if not other:
            return self
        return JsonMessages(self.messages + other.messages)

    def __bool__(self):
        return bool(self.messages)

    def __iter__(self) -> Iterator[JsonMessage]:
        return iter(self.messages)

    def to_json(self):
        return self.messages

    @staticmethod
    def deserialize(json_data: JsonDataType):
        return JsonMessages(messages=[JsonMessage.deserialize(row) for row in json_data])


JSON_MESSAGES_EMPTY = JsonMessages()


@dataclass(frozen=True)
class ValidatedJson:
    """
    ValidatedJson can have a base bit of JSON that's either pure JSON, or other ValidatedJson
    Allowing validation messages to be associated with the parts that caused the problems
    """
    json_data: JsonDataType
    # note that JsonMessgaes are immutable so JSON_MESSAGES_EMPTY can be provided as a default value
    messages: JsonMessages = JSON_MESSAGES_EMPTY
    void: bool = False

    @staticmethod
    def make_void(messages: Optional[JsonMessages] = JSON_MESSAGES_EMPTY):
        return ValidatedJson(
            json_data=None,
            messages=messages,
            void=True
        )

    def to_json(self) -> JsonDataType:
        """
        Unlike most to_json methods, this will lose data as it will return the pure JSON
        stripping away representation of the validation messages.
        Use serialize, deserialize to keep the messages.
        """
        return ValidatedJson.recursive_to_json(self)

    @staticmethod
    def recursive_to_json(data: JsonDataType) -> JsonDataType:
        if isinstance(data, ValidatedJson):
            return ValidatedJson.recursive_to_json(data.json_data)
        elif isinstance(data, dict):
            pure_dict = dict()
            for key, value in data.items():
                if isinstance(value, ValidatedJson) and value.void:
                    pass
                else:
                    pure_dict[key] = ValidatedJson.recursive_to_json(value)
            return pure_dict
        elif isinstance(data, list):
            return [ValidatedJson.recursive_to_json(item) for item in data]
        else:
            return data

    @staticmethod
    def _serialize(obj) -> JsonDataType:
        if isinstance(obj, ValidatedJson):
            if obj.messages or obj.void:
                return {
                    "*wrapper$": "VJ",
                    "messages": obj.messages,
                    "wrap": ValidatedJson._serialize(obj.json_data),
                    "void": obj.void
                }
            else:
                return ValidatedJson._serialize(obj.json_data)
        elif isinstance(obj, list):
            return [ValidatedJson._serialize(entry) for entry in obj]
        elif isinstance(obj, dict):
            return {key: ValidatedJson._serialize(value) for (key, value) in obj.items()}
        else:
            # primitive
            return obj

    def serialize(self) -> JsonDataType:
        return json.loads(json.dumps(ValidatedJson._serialize(self)))

    @staticmethod
    def _deserialize(json_thing: JsonDataType) -> Union[JsonDataType, 'ValidatedJson']:
        if isinstance(json_thing, dict):
            if json_thing.get("*wrapper$") == "VJ":
                return ValidatedJson(
                    json_data=ValidatedJson._deserialize(json_thing.get('wrap')),
                    messages=JsonMessages.deserialize(json_thing.get('messages')),
                    void=json_thing.get('void') == True
                )
            else:
                return {key: ValidatedJson._deserialize(value) for (key, value) in json_thing.items()}
        elif isinstance(json_thing, list):
            return [ValidatedJson._deserialize(entry) for entry in json_thing]
        else:
            # primitive
            return json_thing

    @staticmethod
    def deserialize(json_thing: JsonDataType) -> 'ValidatedJson':
        d_ser = ValidatedJson._deserialize(json_thing)
        if isinstance(d_ser, ValidatedJson):
            return d_ser
        else:
            return ValidatedJson(d_ser)

    def pure_json(self) -> JsonDataType:
        """
        returns dicts, lists, str, int, etc not combined with
        """
        return json.loads(json.dumps(self.to_json()))

    @staticmethod
    def _traverse_messages(json_data) -> JsonMessages:
        messages = JSON_MESSAGES_EMPTY
        if isinstance(json_data, list):
            for val in json_data:
                messages += ValidatedJson._traverse_messages(val)
        elif isinstance(json_data, dict):
            for val in json_data.values():
                messages += ValidatedJson._traverse_messages(val)
        elif isinstance(json_data, ValidatedJson):
            messages += json_data.messages
            messages += ValidatedJson._traverse_messages(json_data.json_data)
        return messages

    def __setitem__(self, key, value):
        lazy.invalidate(self, 'all_messages')
        self.json_data[key] = value

    def __getitem__(self, item):
        return self.json_data[item]

    @lazy
    def all_messages(self) -> JsonMessages:
        """
        Returns messages associated to this ValidatedJson or any of its data
        Where the "messages" property will only pertain to this level
        """
        return ValidatedJson._traverse_messages(self)

    @property
    def has_errors(self) -> bool:
        return any(message.is_error for message in self.all_messages)

    def __copy__(self):
        # messages don't need deeop copy as they're immutable, data does though
        return ValidatedJson(json_data=copy.deepcopy(self.json_data), messages=self.messages)

    def __bool__(self):
        return bool(self.json_data) or bool(self.messages)
