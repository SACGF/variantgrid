from typing import Optional, List, Union
from django import template
from django.db.models import QuerySet

register = template.Library()

@register.simple_tag()
def count(items, singular: str, plural: Optional[str] = None):
    item_count: int
    if isinstance(items, list):
        item_count = len(items)
    elif isinstance(items, QuerySet):
        item_count = items.count()
    else:
        item_count = items

    if item_count == 1:
        return f'{singular}'
    else:
        if not plural:
            plural = f'{singular}s'
        return f'{item_count} {plural}'