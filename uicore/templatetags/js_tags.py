import collections.abc
import datetime
import json
import re
import urllib
from datetime import date, timedelta
from decimal import Decimal
from html import escape
from typing import Union, Any, Optional
import math
from django import template
from django.db.models import TextChoices
from django.utils.safestring import mark_safe, SafeString

from library.utils import format_significant_digits, JsonDataType, format_diff_text
from snpdb.user_settings_manager import UserSettingsManager
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
        text = text.replace('\\', '\\\\').replace('`', '\\`').replace('</script>', '<\\/script>')
        return mark_safe(text)
    return ''

@register.filter
def js_symbol(value):
    """ Used for making a valid symbol (eg function name out of T2T-CHM13v2.0) """
    # Replace invalid characters with underscores
    name = re.sub(r'[^a-zA-Z0-9_$]', '_', value)
    # Ensure it starts with a valid character
    if not re.match(r'^[a-zA-Z_$]', name):
        name = f'_{name}'
    return name


@register.filter
def limit_length(text, limit=100):
    if limit and text and len(text) > limit:
        return text[0:(limit-3)] + '...'
    return text


@register.inclusion_tag("uicore/tags/limit_length_with_tooltip.html")
def limit_length_with_tooltip(text, limit=100):
    tooltip = ""
    if limit and len(text) > limit:
        tooltip = text
        text = text[0:(limit-3)] + '...'
    return {"tooltip": tooltip, "text": text}


@register.filter(is_safe=True)
def format_value_show_invisible(val):
    if isinstance(val, str):
        return format_diff_text(val)
    return format_value(val)


@register.filter(is_safe=True)
def format_value(val, limit=0):
    if val is None:
        return mark_safe('<span class="no-value">None</span>')
    if val == "":
        return mark_safe('<span class="none">""</span>')
    elif isinstance(val, (datetime.datetime, datetime.date)):
        return f"{val:%Y-%m-%d}"
    if isinstance(val, (dict, list)):
        return mark_safe(f'<span class="json">{escape(json.dumps(val))}</span>')
    if isinstance(val, float):
        val = format(Decimal(str(val)).normalize(), 'f')
        return mark_safe(f'<span class="number">{val}</span>')
    if isinstance(val, int):
        if val == 0:
            return mark_safe('<span class="no-value">0</span>')
        return mark_safe(f'<span class="number">{val}</span>')

    if isinstance(val, SafeString):
        return val

    val = str(val)
    if limit:
        val = limit_length(val, limit)

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
    if val is None or len(str(val).strip()) == 0:
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
        # if we're formatting ValidatedJson and the first element has messages, that provides formatting we don't need code-block
        # Also if we already have a css_class (like card-body) we don't need code-block
        if not data or not isinstance(data, dict) or not ('*wrapper$' in data and data.get('messages')):
            css_class = "code-block"
    return {"data": data, "css_class": css_class}


@register.inclusion_tag("uicore/tags/code_block_xml.html")
def code_xml(data: str, css_class: Optional[str] = "code-block"):
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
        # detect if the regex looks like it's just trying to match a file extnesion
        extension = match.group('extension')
    return {"error": error, "extension": extension, "pattern": pattern}


@register.inclusion_tag("uicore/tags/code_shell.html")
def code_shell(data: str):
    # doesn't really do anything of note currently, but gives us the opportunity in future
    return {"text": data}


@register.inclusion_tag("uicore/tags/timedelta.html", name="timedelta")
def timedelta_tag(time: timedelta, show_micro=False):
    remainder = time.seconds
    hours = math.floor(remainder / (60*60))
    if hours:
        remainder -= (hours * 60*60)
    minutes = math.floor(remainder / 60)
    if minutes:
        remainder -= (minutes * 60)
    seconds = f"{remainder:02d}"
    micro = None
    if show_micro:
        micro = time.microseconds

    return {
        "days": time.days,
        "hours": hours,
        "minutes": minutes,
        "seconds": seconds,
        "micro": micro
    }


@register.inclusion_tag("uicore/tags/timestamp.html")
def timestamp(timestamp, time_ago: bool = False, date_only: bool = False, show_seconds: bool = False, show_micro = False, text_only: bool = False, tooltip: str = ""):
    css_classes = []
    if time_ago:
        css_classes.append('time-ago')
    if show_micro:
        css_classes.append('micro')
    elif show_seconds:
        css_classes.append('seconds')
    elif date_only:
        css_classes.append('date-only')

    date_value = None
    if timestamp:
        if not isinstance(timestamp, (int, float)):
            if not hasattr(timestamp, 'timestamp'):
                if isinstance(timestamp, date):
                    date_value = timestamp
                    timestamp = datetime.datetime(year=timestamp.year, month=timestamp.month, day=timestamp.day, tzinfo=UserSettingsManager.get_user_timezone())
                else:
                    raise ValueError(f"Unsure how to convert {timestamp} to timestamp")
            timestamp = timestamp.timestamp()
        return {
            "tooltip": tooltip,
            "datetime": datetime.datetime.fromtimestamp(timestamp),
            "date": date_value,
            "timestamp": timestamp,
            "css_class": " ".join(css_classes),
            "text_only": text_only
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
    elif isinstance(value, TextChoices):
        return value.label
    return value


# Filters to make up for how purposefully crippled Django Templates are


@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)


@register.filter
def times(value):
    return range(0, value)
