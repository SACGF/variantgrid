"""
Template tag for naming the classify queue vocabulary in help text, so a page says the tags this
deployment actually has rather than a name baked into the copy.

Entry point: {% classify_queue_tag_names %} - @see snpdb/models/models.py:Tag.classify_queue_qs
"""
from django import template

from library.utils.text_utils import join_with_commas_and_ampersand
from snpdb.models import Tag

register = template.Library()


@register.simple_tag
def classify_queue_tag_names() -> str:
    """ The queue tag names, ready to drop into a sentence ('ToDo or SomaticToDo'). Empty when nothing is
        flagged as a queue tag, so help text can leave the whole tagging clause out """
    tag_ids = list(Tag.classify_queue_qs().order_by("pk").values_list("pk", flat=True))
    return join_with_commas_and_ampersand(tag_ids, final_sep="or")
