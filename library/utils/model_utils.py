from dataclasses import dataclass
from typing import Type, Callable, Optional, TypeVar, Generic, Union
from django.core.exceptions import FieldDoesNotExist
from django.db import models
from django.db.models import Model
import json


class ModelUtilsMixin:
    """
    Allows you to deal with an instance of a model or the key of the model
    """

    @classmethod
    def get(cls: Type[Model], value):
        if value is None:
            return None
        if isinstance(value, cls):
            return value
        if isinstance(value, (str, int)):
            return cls.objects.get(pk=value)
        raise ValueError(f'Expected {cls.__name__} or str or int, got {value}')


class ArrayLength(models.Func):
    """
    Can annotate with array length now e.g.
    MyModel.objects.all().annotate(field_len=ArrayLength('field')).order_by('field_len')
    """
    function = 'CARDINALITY'


class AsciiValue(models.Func):
    function = 'ASCII'


def model_has_field(model: Type[Model], field_name: str) -> bool:
    try:
        if field_name.endswith('_id'):
            field = model._meta.get_field(field_name.strip('_id'))
            if field.is_relation:
                return True
        else:
            _ = model._meta.get_field(field_name)
            return True
    except FieldDoesNotExist:
        pass

    return False


T = TypeVar("T")


@dataclass
class AuditSingleChange(Generic[T]):
    value: T
    log_entry: 'LogEntry'

    @property
    def user(self):
        return self.log_entry.actor

    @property
    def timestamp(self):
        return self.log_entry.timestamp


class AuditUtils:

    @staticmethod
    def last_change_for(model_instance: Model, field: str, is_json: bool = False, parser: Optional[Callable[[Union[str, dict]], T]] = None) -> AuditSingleChange[T]:
        from auditlog.models import LogEntry
        if log_entry := (LogEntry.objects.get_for_object(model_instance)
                .filter(**{f"changes__{field}__isnull": False})
                .exclude(**{f"changes__{field}__1": "None"})  # the changes are stored very stringified, to the point where None is saved as "None"
                .order_by('-timestamp').first()):
            value = log_entry.changes.get(field)[1]
            if isinstance(value, str) and is_json:
                try:
                    value = json.loads(value)
                    if isinstance(value, str):
                        # not sure what's going on with double encoding
                        value = json.loads(value)
                except:
                    pass
            if parser:
                try:
                    value = parser(value)
                except Exception as ex:
                    raise ex
                    # raise ValueError(f"Couldn't parse {field} \"{value_str}\"")
            return AuditSingleChange(value, log_entry)
        return None, None
