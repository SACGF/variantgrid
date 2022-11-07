from typing import Type
from django.core.exceptions import FieldDoesNotExist
from django.db import models
from django.db.models import Model


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
