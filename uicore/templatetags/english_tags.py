from typing import Optional, List, Union

from django import template
from django.db.models import QuerySet

from library.utils import pretty_label

register = template.Library()


def count_items(items: Union[List, QuerySet, int]) -> int:
    item_count: int
    if hasattr(items, '__len__'):
        item_count = len(items)
    elif isinstance(items, QuerySet):
        item_count = items.count()
    else:
        item_count = items

    return item_count


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
    item_count = count_items(items)

    if item_count == 1:
        return f'{singular}'
    else:
        if plural is None:
            plural = f'{singular}s'
        return f'{item_count} {plural}'


@register.simple_tag()
def plural(items: Union[List, QuerySet, int], singular: str = "", plural: str = "s"):
    if count_items(items) == 1:
        return singular
    return plural


@register.filter()
def code_to_english(text: str):
    return pretty_label(text)
