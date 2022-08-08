import operator
from functools import reduce
from typing import Optional, Any

from django.db.models import Q


class QueryJsonFilter:
    """
    Allows you to generate a filter based on JSON, limited abilities

    "key": "value" (does exact)
    "key": True (makes sure key is not null)
    "key": [value1, value2, valueN] (makes sure key is not null and is one of the provided values)
    "key": [value1, value2, valueN, null] (makes sure key is one of the provided values including null)
    "key": ["!", value1, value2, valueN] (makes sure key is not one of the provided values, but can be null)
    "key": ["!", value1, value2, valueN, null] (makes sure key is not any of the provided values and can not be null)
    "key": {operation: value} (lets you run other queries, e.g. operation can be "lt", "lte", "iexact", "contains"}
    "or": [queries] (ors items together, if multiple queries provided without "or", and is used by default)
    "and": [queries] (ands items together)
    "not": {queries} (negate queries)
    """

    def __init__(self, prefix: str, postfix: str):
        self.prefix = prefix
        self.postfix = postfix

    @staticmethod
    def classification_value_filter() -> 'QueryJsonFilter':
        return QueryJsonFilter(prefix='published_evidence__', postfix='__value')

    def q(self, key: str, value: Any, operation: str = 'exact') -> Q:
        return Q(**{f"{self.prefix}{key}{self.postfix}__{operation}": value})

    def convert_to_q(self, blob, op=operator.__and__) -> Q:
        if isinstance(blob, list):
            return reduce(op, (self.convert_to_q(entry) for entry in blob))
        elif isinstance(blob, dict):
            return reduce(op, (self.convert_to_q_w_key(key, value) for key, value in blob.items()))
        else:
            raise ValueError(f"Expected list or dictionary, found {blob}")

    def convert_to_q_w_key(self, key: str, value) -> Q:
        # special hardcoded segment, should migrate
        if key == 'somatic':
            if value is False:
                # old code also allowed value to be isnull
                return ~self.q('allele_origin', ['somatic'], 'in') | self.q('allele_origin', True, 'isnull')
            else:
                # so ugly that this is the only way I know to make a Q that wont restrict record counts
                raise ~Q(pk=None)
        else:
            if key == 'not':
                return ~self.convert_to_q(value)
            elif key == 'or':
                return self.convert_to_q(value, operator.__or__)
            elif key == 'and':
                return self.convert_to_q(value, operator.__and__)
            else:
                # key is assumed to be regular value key
                if isinstance(value, list):
                    handle_none = False
                    is_not = False
                    if value[0] == '!':
                        # treat this as not_in
                        value = value[1:]
                        is_not = True
                    handle_none = False
                    if None in value:
                        value.remove(None)
                        handle_none = True
                    q = self.q(key, value, 'in')
                    if is_not:
                        q = ~q
                        if handle_none:
                            q = q & self.q(key, False, 'isnull')
                        else:
                            # asking for allele origin to not be somatic, still exclude nulls
                            # so we want ["!", "somatic"] to allow null and all non somatic values
                            q = q | self.q(key, True, 'isnull')
                    elif handle_none:
                        q = q | self.q(key, True, 'isnull')

                    return q
                elif value is None:
                    return self.q(key, True, 'isnull')
                elif isinstance(value, dict):
                    if len(value) != 1:
                        raise ValueError(f"{key} can only have a single entry dict for value, but found {value}")
                    for sub_key, sub_value in value.items():
                        # should only ever be one loop here
                        return self.q(key, value=sub_value, operation=sub_key)
                else:
                    return self.q(key, value)
