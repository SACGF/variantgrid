"""
From https://djangosnippets.org/snippets/2696/
"""

from django import forms
from django.forms.widgets import TextInput

class NumberInput(TextInput):
    input_type = 'number'


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
