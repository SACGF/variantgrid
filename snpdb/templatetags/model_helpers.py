"""

By Fydo from http://stackoverflow.com/a/26614950

"""
from django import template
from django.utils.html import escape

from library.django_utils import get_model_fields_and_formatted_values_tuples_list

register = template.Library()


def _display_rows(model):
    """ (label, css data type, escaped value) per model field """
    for name, field in get_model_fields_and_formatted_values_tuples_list(model):
        data_type = field.__class__.__name__
        name = name.replace('_', ' ')

        if str(field).isdigit():
            data_type = 'num'
            value = format(int(field), ',')
        else:
            value = escape(field)
        yield name, data_type, value


@register.filter()
def as_table(model):
    # Callers apply |safe at the template, so this stays a plain str
    return "".join(f'<tr id="{name}-row"><th class="name {data_type}">{name}</th>'
                   f'<td class="field {data_type}">{value}</td></tr>'
                   for name, data_type, value in _display_rows(model))


@register.filter()
def display_rows(model):
    """ For a {% labelled %} row per model field """
    return list(_display_rows(model))
