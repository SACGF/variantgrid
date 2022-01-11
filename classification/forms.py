from django import forms

from classification.models import EvidenceKey
from library.django_utils.autocomplete_utils import ModelSelect2


class EvidenceKeyForm(forms.Form):
    evidence_key = forms.ModelChoiceField(queryset=EvidenceKey.objects.all(),
                                          required=False,
                                          widget=ModelSelect2(url='evidence_key_autocomplete',
                                                              attrs={'data-placeholder': 'EKey...',
                                                                     'class': 'custom-key'}))
    value = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'value...', 'class': 'custom-value'}),
                            required=False)
