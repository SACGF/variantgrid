"""
From https://djangosnippets.org/snippets/2696/
"""

import itertools

from django import forms
from django.forms.widgets import TextInput, Widget
from django.utils.html import format_html
from django.utils.safestring import mark_safe


class NumberInput(TextInput):
    input_type = 'number'


class StarsWidget(Widget):
    """ Interactive FontAwesome star-rating input (replacement for the django_starfield 'Stars' widget).

        Renders hidden radio inputs styled as stars via css/stars_input.css. Filled/hover behaviour is
        pure CSS (the `input:checked ~ label` sibling selector), so the stars are rendered highest-first.
        Matches the read-only `clinvar_stars` template tag styling (fa-solid star + text-success).

        A zero radio is rendered last, labelled "0" - it sits to the left of the 1 star (the row is
        laid out right-to-left) so "no minimum" is picked the same way as any other value, and always
        submits a value - an omitted radio group fails validation on a required IntegerField. """

    class Media:
        css = {"all": ["css/stars_input.css"]}

    def __init__(self, attrs=None, stars=4):
        super().__init__(attrs)
        self.stars = stars

    def render(self, name, value, attrs=None, renderer=None):
        container_id = f"stars-{name}"
        parts = []
        # Highest first so the CSS sibling selector can fill the selected star and all lower ones.
        for i in range(self.stars, 0, -1):
            input_id = f"{name}-{i}"
            checked = mark_safe(" checked") if str(value) == str(i) else ""
            title = f"{i} star" + ("" if i == 1 else "s")
            parts.append(format_html(
                '<input type="radio" name="{name}" id="{id}" value="{i}"{checked}>'
                '<label for="{id}" title="{title}"><i class="fa-regular fa-star"></i></label>',
                name=name, id=input_id, i=i, checked=checked, title=title))

        zero_id = f"{name}-0"
        zero_checked = mark_safe(" checked") if str(value) not in [str(i) for i in range(1, self.stars + 1)] else ""
        parts.append(format_html(
            '<input type="radio" name="{name}" id="{id}" value="0"{checked}>'
            '<label class="stars-none" for="{id}" title="No minimum">0</label>',
            name=name, id=zero_id, checked=zero_checked))
        return format_html('<span class="stars-input-container">'
                           '<span class="stars-input" id="{container_id}">{stars}</span>'
                           '</span>',
                           container_id=container_id, stars=mark_safe("".join(parts)))


class ReadOnlyDisplayWidget(Widget):
    """ Renders the value as plain text instead of a form control - a disabled Select still
        renders every option, so use this for choice fields with a large pool (eg users) """

    def __init__(self, display, attrs=None):
        super().__init__(attrs)
        self.display = display

    def render(self, name, value, attrs=None, renderer=None):
        return format_html('<span class="form-control-plaintext">{}</span>', self.display)


class ROFormMixin(forms.BaseForm):
    """ Default is make all fields ready only. Use "read_only" tuple for a subset of fields.
        Fields in "read_only_display" are also disabled, and render as plain text """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        read_only_display_fields = getattr(self.Meta, "read_only_display", ())
        read_only_fields = getattr(self.Meta, "read_only", None)
        if read_only_fields is None and not read_only_display_fields:
            read_only_fields = self.fields  # all

        for field in itertools.chain(read_only_fields or (), read_only_display_fields):
            self.fields[field].disabled = True

        instance = getattr(self, "instance", None)
        for field_name in read_only_display_fields:
            value = getattr(instance, field_name, None)
            self.fields[field_name].widget = ReadOnlyDisplayWidget(str(value) if value is not None else "")
