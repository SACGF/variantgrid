import re
import uuid
from typing import Any

from django.template import Library

register = Library()

just_zeros = re.compile('^0[.]0+$')
"""
@register.inclusion_tag("snpdb/tags/labelled.html")
def labelled(id: str = None, label: str = None, value: Any = None, css_class: str = None, span_class: str = None, title: str = None):
    if not id:
        id = str(uuid.uuid4())
    if not label:
        label = id.replace('_', ' ')
    if isinstance(value, str):
        value = value.strip()
        if value == '' or value == '-':
            value = None
        elif just_zeros.match(value):
            span_class = 'zero-value'
    return {"id": id, "label": label, "value": value, "css_class": css_class, "span_class": span_class, "title": title}
"""
@register.inclusion_tag("snpdb/tags/labelled.html")
def labelled_old(id: str = None, label: str = None, value: Any = None, css_class: str = None, span_class: str = None, title: str = None):
    if not id:
        id = str(uuid.uuid4())
    if not label:
        label = id.replace('_', ' ')
    if isinstance(value, str):
        value = value.strip()
        if value == '' or value == '-':
            value = None
        elif just_zeros.match(value):
            span_class = 'zero-value'
    return {"id": id, "label": label, "value": value, "css_class": css_class, "span_class": span_class, "title": title}
