"""

By Fydo from http://stackoverflow.com/a/26614950

"""
from django import template
from django.utils.html import escape
from django.utils.safestring import mark_safe

from library.django_utils import get_model_fields, get_model_fields_and_formatted_values_tuples_list

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
def as_p(model):
    return mark_safe("".join(f'<p><label for="{name}-value">{name}</label>'
                             f'<span class="field {data_type}">{value}</span></p>'
                             for name, data_type, value in _display_rows(model)))

@register.filter()
def qs_as_htable(qs):
    ret = ""
    data = list(qs)
    if data:
        model = data[0]
        ret = "<table>"
        ret += "<tr>"

        field_names = get_model_fields(model)
        for name in field_names:
            ret += '<th class="name">' + name + '</th>'

        for model in data:
            ret += "<tr>"
            for name in field_names:
                try:
                    get_display_func = getattr(model, f"get_{name}_display")
                    field = get_display_func()
                except:
                    field = str(getattr(model, name))

                if not field:
                    field = ''
                ret += '<td class="name">' + field + '</td>'

            ret += "</tr>"
        ret += "</table>"
    return ret
