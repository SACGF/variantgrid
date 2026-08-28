""" Column filter rules - the vocabulary shared by everything that filters a grid by column value.

A rule set is ``{"groupOp": "AND"|"OR", "rules": [{"field": ..., "op": ..., "data": ...}]}``. It is
what the DataTables filter builder sends up, what FilterNode persists as FilterNodeItem rows, and
what filter_rules_to_q turns into a Django Q - so all three stay in step off FILTER_OPERATIONS.
"""
import json
import operator
from functools import reduce
from typing import Optional

from django.db.models.query_utils import Q
from django.utils.encoding import smart_str

from library.utils import JsonObjType

# 'takes_data' is False for the operations that filter on their own
FILTER_OPERATIONS = [
    ("eq", "equal", True),
    ("ne", "not equal", True),
    ("bw", "begins with", True),
    ("bn", "does not begin with", True),
    ("ew", "ends with", True),
    ("en", "does not end with", True),
    ("cn", "contains", True),
    ("nc", "does not contain", True),
    ("nu", "is null", False),
    ("nn", "is not null", False),
    ("in", "is in", True),
    ("ni", "is not in", True),
    ("gt", "greater than", True),
    ("ge", "greater than or equal to", True),
    ("lt", "less than", True),
    ("le", "less than or equal to", True),
]

FILTER_OPERATION_LABELS = {op: label for op, label, _takes_data in FILTER_OPERATIONS}

# operation -> (django lookup, exclude)
_CASE_SENSITIVE_LOOKUPS = {
    'ne': ('%(field)s__exact', True),
    'bn': ('%(field)s__startswith', True),
    'en': ('%(field)s__endswith', True),
    'nc': ('%(field)s__contains', True),
    'ni': ('%(field)s__in', True),
    'in': ('%(field)s__in', False),
    'eq': ('%(field)s__exact', False),
    'bw': ('%(field)s__startswith', False),
    'gt': ('%(field)s__gt', False),
    'ge': ('%(field)s__gte', False),
    'lt': ('%(field)s__lt', False),
    'le': ('%(field)s__lte', False),
    'ew': ('%(field)s__endswith', False),
    'cn': ('%(field)s__contains', False),
    'nu': ('%(field)s__isnull', True),
    'nn': ('%(field)s__isnull', False),
}

_CASE_INSENSITIVE_LOOKUPS = _CASE_SENSITIVE_LOOKUPS | {
    'ne': ('%(field)s__iexact', True),
    'eq': ('%(field)s__iexact', False),
    'bn': ('%(field)s__istartswith', True),
    'bw': ('%(field)s__istartswith', False),
    'en': ('%(field)s__iendswith', True),
    'ew': ('%(field)s__iendswith', False),
    'nc': ('%(field)s__icontains', True),
    'cn': ('%(field)s__icontains', False),
}


def format_operation(op: str) -> str:
    return FILTER_OPERATION_LABELS[op]


def filter_operations() -> list[JsonObjType]:
    """ The operations the client's filter builder offers """
    return [{"op": op, "label": label, "takesData": takes_data} for op, label, takes_data in FILTER_OPERATIONS]


def filter_rules_from_params(params) -> Optional[JsonObjType]:
    """ Rules off a request's GET/POST params, or None where it carries none """
    if rules_json := params.get("filters"):
        try:
            return json.loads(rules_json)
        except ValueError:
            return None
    return None


def filter_rules_to_q(rules: Optional[JsonObjType], ignore_case: bool = True) -> Optional[Q]:
    # TODO: Add more support for RelatedFields (searching and displaying)
    # FIXME: Validate data types are correct for field being searched.
    if not rules:
        return None

    lookups = _CASE_INSENSITIVE_LOOKUPS if ignore_case else _CASE_SENSITIVE_LOOKUPS
    q_filters = []
    for rule in rules['rules']:
        op, field, data = rule['op'], rule['field'], rule['data']

        filter_fmt, exclude = lookups[op]
        filter_str = smart_str(filter_fmt % {'field': field})
        if filter_fmt.endswith('__in'):
            filter_kwargs = {filter_str: data.split(',')}
        elif filter_fmt.endswith('__isnull'):
            # FilterNode was slow - pass exclude as arg and don't invert Q
            # Prev generated code like:
            # (NOT (AND: ('variantannotation__dbsnp_rs_id__isnull', True))
            # which generated a full table scan on variant (>100M+ rows...)
            filter_kwargs = {filter_str: exclude}
            exclude = False
        else:
            filter_kwargs = {filter_str: smart_str(data)}

        q = Q(**filter_kwargs)
        if exclude:
            q = ~q
        q_filters.append(q)

    if not q_filters:  # Filter w/no rules
        return None
    if rules['groupOp'].upper() == 'OR':
        return reduce(operator.or_, q_filters)
    return reduce(operator.and_, q_filters)
