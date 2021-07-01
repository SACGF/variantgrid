import urllib
from decimal import Decimal
from html import escape
from typing import Union, Any, Optional

from django import template
from django.utils.safestring import mark_safe
import json
import re

# FIXME, move this out of classifications and into snpdb
from classification.views.classification_datatables import DatatableConfig
from library.utils import format_significant_digits

register = template.Library()


def jsonify_for_js(json_me, pretty=False) -> Union[str, Any]:
    if isinstance(json_me, str):
        return json_me
    elif isinstance(json_me, bool):
        if json_me:
            return mark_safe('true')
        else:
            return mark_safe('false')
    elif isinstance(json_me, int) or isinstance(json_me, float):
        return json_me
    indent = 0 if not pretty else 4
    text = json.dumps(json_me, indent=indent)
    if pretty:
        # this stops arrays of arrays taking up too much vertical space
        text = re.compile(r'],\s*\[', re.MULTILINE).sub('],[', text)
    text = text.replace('</script>', '<\\/script>')
    return mark_safe(text)


@register.filter
def jsonify(json_me) -> Union[str, Any]:
    return jsonify_for_js(json_me)


@register.filter
def jsonify_pretty(json_me) -> Union[str, float]:
    return jsonify_for_js(json_me, pretty=True)

@register.filter
def query_unquote(query_string):
    return urllib.parse.unquote(query_string)


@register.filter
def jsstring(text):
    if text:
        text = text.replace('\\', '\\\\').replace('`','\`').replace('</script>','<\\/script>')
        return mark_safe(text)
    else:
        return ''


@register.filter
def limit_length(text, limit=100):
    if text and len(text) > limit:
        return text[0:(limit-3)] + '...'
    return text


@register.filter(is_safe=True)
def format_value(val):
    if val is None:
        return mark_safe('<span class="no-value">None</span>')
    if val == "":
        return mark_safe('<span class="none">""</span>')
    if isinstance(val, dict) or isinstance(val, list):
        return mark_safe(f'<span class="json">{escape(json.dumps(val))}</span>')
    if isinstance(val, float):
        val = format(Decimal(str(val)).normalize(), 'f')
        return mark_safe(f'<span class="number">{val}</span>')
    elif isinstance(val, int):
        return mark_safe(f'<span class="number">{val}</span>')
    else:
        val = escape(val)
        val = val.replace("\n", "<br/>")
        return mark_safe(val)


@register.filter()
def format_computer_text(val):
    if val is None:
        return ''
    return val.replace('&', ' & ').replace('_', ' ')


@register.filter(is_safe=True)
def format_unit_as_percent(val: Optional[float]):
    if val is None:
        return ''
    return f"{format_significant_digits(val*100)} %"


@register.filter()
def dash_if_empty(val):
    if val is None or len(val.strip()) == 0:
        return mark_safe('<span class="no-value">-</span>')
    return val


@register.inclusion_tag("classification/tags/timestamp.html")
def timestamp(timestamp, time_ago: bool = False):
    css_class = 'time-ago' if time_ago else ''
    if timestamp:
        if not isinstance(timestamp, int) and not isinstance(timestamp, float):
            timestamp = timestamp.timestamp()
        return {
            "timestamp": timestamp,
            "css_class": css_class
        }
    else:
        return {"css_class": "empty"}


@register.filter()
def duration(td):
    total_seconds = int(td.total_seconds())

    days = total_seconds // 86400
    remaining_hours = total_seconds % 86400
    remaining_minutes = remaining_hours % 3600
    hours = remaining_hours // 3600
    minutes = remaining_minutes // 60
    seconds = remaining_minutes % 60

    days_str = f'{days}d ' if days else ''
    hours_str = f'{hours}h ' if hours else ''
    minutes_str = f'{minutes}m ' if minutes else ''
    seconds_str = f'{seconds}s' if seconds and not hours_str else ''

    return f'{days_str}{hours_str}{minutes_str}{seconds_str}'


@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)


@register.filter
def format_preference(value):
    if value is True:
        return 'Yes'
    elif value is False:
        return 'No'
    else:
        return value


@register.filter
def times(value):
    return range(0, value)
