"""
From https://djangosnippets.org/snippets/2696/
"""

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
        Matches the read-only `clinvar_stars` template tag styling (fa-solid star + text-success). """

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
        return format_html('<div class="stars-input" id="{}">{}</div>',
                           container_id, mark_safe("".join(parts)))


class ROFormMixin(forms.BaseForm):
    """ Default is make all fields ready only. Use "read_only" tuple for a subset of fields  """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if hasattr(self.Meta, "read_only"):
            read_only_fields = self.Meta.read_only
        else:
            read_only_fields = self.fields  # all

        for field in read_only_fields:
            self.fields[field].disabled = True
