"""
Column filter rules - the vocabulary shared by the grids' filter builder and FilterNode.

A rule set looks like {'groupOp': 'AND'|'OR', 'rules': [{'field': ..., 'op': ..., 'data': ...}]}:
the JSON the DataTables filter builder emits, what FilterNodeItem rows are written from, and what
the grid data endpoint reads off the 'filters' param. rules_to_q turns it into the Q object both
paths filter on, so a grid column filter and a FilterNode on the same rules produce the same
queryset.
"""
import json
import logging
import operator
from functools import reduce
from typing import Optional

from django.db.models.query_utils import Q
from django.utils.encoding import smart_str

from library.utils import JsonObjType

logger = logging.getLogger(__name__)

# rules_to_q below turns these into Django lookups, FilterNodeItem persists them, and the DataTables
# filter builder offers them - so all three stay in step off this one list.
# 'takes_data' is False for the operations that filter on their own.
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

# op: (django_lookup, use_exclude). Text comparisons are case insensitive
_FILTER_MAP = {
    'ne': ('%(field)s__iexact', True),
    'bn': ('%(field)s__istartswith', True),
    'en': ('%(field)s__iendswith', True),
    'nc': ('%(field)s__icontains', True),
    'ni': ('%(field)s__in', True),
    'in': ('%(field)s__in', False),
    'eq': ('%(field)s__iexact', False),
    'bw': ('%(field)s__istartswith', False),
    'gt': ('%(field)s__gt', False),
    'ge': ('%(field)s__gte', False),
    'lt': ('%(field)s__lt', False),
    'le': ('%(field)s__lte', False),
    'ew': ('%(field)s__iendswith', False),
    'cn': ('%(field)s__icontains', False),
    'nu': ('%(field)s__isnull', True),
    'nn': ('%(field)s__isnull', False),
}


def format_operation(op: str) -> str:
    return FILTER_OPERATION_LABELS[op]


def filter_operations_json() -> list[JsonObjType]:
    """ The operations the filter builder offers """
    return [{"op": op, "label": label, "takesData": takes_data} for op, label, takes_data in FILTER_OPERATIONS]


def parse_filters(param: Optional[str]) -> Optional[dict]:
    """ The 'filters' request param as a rule set, or None if absent/unparseable """
    if not param:
        return None
    try:
        filters = json.loads(param)
    except ValueError:
        logger.warning("Could not parse filter rules: %s", param)
        return None
    if isinstance(filters, dict) and filters.get("rules"):
        return filters
    return None


def rules_to_q(json_filters: dict) -> Optional[Q]:
    q_filters = []
    for rule in json_filters['rules']:
        op, field, data = rule['op'], rule['field'], rule['data']

        filter_fmt, exclude = _FILTER_MAP[op]
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

    if not q_filters:  # Check for search w/no rules
        return None
    if json_filters['groupOp'].upper() == 'OR':
        return reduce(operator.or_, q_filters)
    return reduce(operator.and_, q_filters)
