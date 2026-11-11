"""
Model-level helpers: ModelUtilsMixin (accept an instance or its pk), the ArrayLength / AsciiValue SQL
functions, model_has_field, refresh_for_update (re-read a row under SELECT FOR UPDATE),
get_model_content_type_dict, and AuditUtils.last_change_for / last_changes_for (a field's latest value from
django-auditlog, for one record or many).
"""
import datetime
import json
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Callable, Generic, Optional, TypeVar, Union

from auditlog.models import LogEntry
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import FieldDoesNotExist
from django.db import models
from django.db.models import Model


class ModelUtilsMixin:
    """
    Allows you to deal with an instance of a model or the key of the model
    """

    @classmethod
    def get(cls: type[Model], value):
        if value is None:
            return None
        if isinstance(value, cls):
            return value
        if isinstance(value, (str, int)):
            return cls.objects.get(pk=value)
        raise ValueError(f'Expected {cls.__name__} or str or int, got {value}')


class ArrayLength(models.Func):
    """
    Annotate with array length e.g.
    MyModel.objects.all().annotate(field_len=ArrayLength('field')).order_by('field_len')
    """
    function = 'CARDINALITY'


class AsciiValue(models.Func):
    function = 'ASCII'


def model_has_field(model: type[Model], field_name: str) -> bool:
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


ModelT = TypeVar("ModelT", bound=Model)


def refresh_for_update(obj: ModelT) -> ModelT:
    return type(obj).objects.filter(pk=obj.pk).select_for_update().get()


def get_model_content_type_dict(model):
    content_type = ContentType.objects.get_for_model(model)
    return {
        'app_label': content_type.app_label,
        'model': content_type.model
    }


T = TypeVar("T")


@dataclass
class AuditSingleChange(Generic[T]):
    value: T
    log_entry: LogEntry

    @property
    def user(self):
        return self.log_entry.actor

    @property
    def timestamp(self) -> datetime:
        return self.log_entry.timestamp


class AuditUtils:

    @staticmethod
    def _single_change(log_entry: LogEntry, field: str, is_json: bool, parser: Optional[Callable[[Union[str, dict]], T]]) -> AuditSingleChange[T]:
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
            value = parser(value)
        return AuditSingleChange(value, log_entry)

    @staticmethod
    def last_change_for(model_instance: Model, field: str, is_json: bool = False, parser: Optional[Callable[[Union[str, dict]], T]] = None) -> AuditSingleChange[T]:
        # order_by = '-timestamp'
        # because we've done some wonky things with timestamp, use the more objective count index for getting the most recent comment
        order_by = '-pk'

        if log_entry := LogEntry.objects.get_for_object(model_instance).filter(**{f"changes__{field}__isnull": False}).order_by(order_by).first():
            return AuditUtils._single_change(log_entry, field, is_json, parser)
        return None

    @staticmethod
    def last_changes_for(model: type[Model], pks: Iterable[int], field: str, is_json: bool = False, parser: Optional[Callable[[Union[str, dict]], T]] = None) -> dict[int, AuditSingleChange[T]]:
        """ last_change_for of many records of one model (integer pks) in one query, keyed by pk - a record without a change to field is absent """
        log_entries = LogEntry.objects.filter(
            content_type=ContentType.objects.get_for_model(model),
            object_id__in=pks,
            **{f"changes__{field}__isnull": False}
        ).order_by('object_id', '-pk').distinct('object_id').select_related('actor')
        return {log_entry.object_id: AuditUtils._single_change(log_entry, field, is_json, parser) for log_entry in log_entries}
