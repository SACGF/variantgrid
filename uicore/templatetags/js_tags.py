import collections.abc
import json
import re
import urllib
from decimal import Decimal
from html import escape
from typing import Union, Any, Optional

from django import template
from django.utils.safestring import mark_safe

from library.utils import format_significant_digits
from uicore.json.json_types import JsonDataType
from uicore.json.validated_json import ValidatedJson

register = template.Library()


def jsonify_for_js(json_me, pretty=False) -> Union[str, Any]:
    if isinstance(json_me, str):
        json_me = json_me.replace('"', '\"')
        return mark_safe(f"\"{json_me}\"")
    if isinstance(json_me, bool):
        if json_me:
            return mark_safe('true')
        return mark_safe('false')
    if isinstance(json_me, (int, float)):
        return json_me
    indent = None if not pretty else 4
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
        text = text.replace('\\', '\\\\').replace('`', '\`').replace('</script>', '<\\/script>')
        return mark_safe(text)
    return ''


@register.filter
def limit_length(text, limit=100):
    if limit and text and len(text) > limit:
        return text[0:(limit-3)] + '...'
    return text


@register.filter(is_safe=True)
def format_value(val):
    if val is None:
        return mark_safe('<span class="no-value">None</span>')
    if val == "":
        return mark_safe('<span class="none">""</span>')
    if isinstance(val, (dict, list)):
        return mark_safe(f'<span class="json">{escape(json.dumps(val))}</span>')
    if isinstance(val, float):
        val = format(Decimal(str(val)).normalize(), 'f')
        return mark_safe(f'<span class="number">{val}</span>')
    if isinstance(val, int):
        return mark_safe(f'<span class="number">{val}</span>')

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


@register.inclusion_tag("uicore/tags/code_block_json.html")
def code_json(data: JsonDataType, css_class: Optional[str] = "", dash_if_empty: bool = False):
    if isinstance(data, ValidatedJson):
        data = data.serialize()

    # note that we still want to print data if it's "false" or "0" but not if it's None or an empty dict or list
    if dash_if_empty and data is None or (isinstance(data, collections.abc.Sized) and len(data) == 0):
        return {"blank": True}

    if not css_class:
        # if we're formatting ValidatedJson and the first element has messages, that provides formatting
        # so we don't need code-block
        # Also if we already have a css_class (like card-body) we don't need code-block
        if not data or not isinstance(data, dict) or ('*wrapper$' not in data and not data.get('messages')):
            css_class = "code-block"
    return {"data": data, "css_class": css_class}


@register.inclusion_tag("uicore/tags/code_block_regex.html")
def code_regex(data: str):
    error = None
    extension = None
    pattern = data
    try:
        re.compile(data)
    except ValueError:
        error = f"Invalid pattern: {data}"

    if match := re.compile(r"\.\*\\\.(?P<extension>.+)").match(data):
        extension = match.group('extension')
    return {"error": error, "extension": extension, "pattern": pattern}


@register.inclusion_tag("uicore/tags/code_shell.html")
def code_shell(data: str):
    # doesn't really do anything of note currently, but gives us the opportunity in future
    return {"text": data}


@register.inclusion_tag("uicore/tags/timestamp.html")
def timestamp(timestamp, time_ago: bool = False, show_seconds: bool = False):
    css_classes = list()
    if time_ago:
        css_classes.append('time-ago')
    if show_seconds:
        css_classes.append('seconds')

    if timestamp:
        if not isinstance(timestamp, (int, float)):
            timestamp = timestamp.timestamp()
        return {
            "timestamp": timestamp,
            "css_class": " ".join(css_classes)
        }
    return {}


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
def format_preference(value):
    if value is True:
        return 'Yes'
    if value is False:
        return 'No'
    return value


# Filters to make up for how purposefully crippled Django Templates are


@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)


@register.filter
def times(value):
    return range(0, value)
