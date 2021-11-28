from django import forms
from django.forms.models import modelform_factory, ALL_FIELDS

from library.forms import ROFormMixin
from .models import CuffDiffFile

EXPRESSION_WIDGETS = {'name': forms.TextInput(),
                      'sample_1': forms.TextInput(),
                      'sample_2': forms.TextInput()}

cdf_model = modelform_factory(CuffDiffFile, fields=ALL_FIELDS, widgets=EXPRESSION_WIDGETS)
class ExpressionFileForm(cdf_model, ROFormMixin):
    pass
