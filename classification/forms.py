from dal import autocomplete
from django import forms

from classification.models import EvidenceKey


class EvidenceKeyForm(forms.Form):
    evidence_key = forms.ModelChoiceField(queryset=EvidenceKey.objects.all(),
                                          required=False,
                                          widget=autocomplete.ModelSelect2(url='evidence_key_autocomplete',
                                                                           attrs={'data-placeholder': 'EKey...',
                                                                                  'class': 'custom-key'}))
    value = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'value...', 'class': 'custom-value'}), required=False)
