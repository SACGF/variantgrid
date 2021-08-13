from typing import Optional, List, Union
from django import template
from django.db.models import QuerySet

register = template.Library()


@register.simple_tag()
def count(items: Union[List, QuerySet, int], singular: str, plural: Optional[str] = None):
    """
    Possibly not the best name, but allows you to change the English on the basis of if you're dealing with singular or plural
    so you don't have to do "item(s)"
    :param items: How many items do you have
    :param singular: How would you refer to it if there's exactly 1, e.g. "one car" to say "Look at my one car"
    :param plural: Refer to it if there's any other number, will be prefixed by the count e.g. "cars" to say "Look at my 2 cars"
    If plural is not provided will default to singular with a s on the end.
    """
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
