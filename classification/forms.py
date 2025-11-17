from django import forms

from classification.models import EvidenceKey
from library.django_utils.autocomplete_utils import ModelSelect2
from uicore.utils.form_helpers import form_helper_horizontal


class EvidenceKeyForm(forms.Form):
    evidence_key = forms.ModelChoiceField(queryset=EvidenceKey.objects.all(),
                                          required=False,
                                          widget=ModelSelect2(url='evidence_key_autocomplete',
                                                              attrs={'data-placeholder': 'EKey...',
                                                                     'class': 'custom-key'}))
    value = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'value...', 'class': 'custom-value'}),
                            required=False)


class ClassificationAlleleOriginForm(forms.Form):
    ALLELE_ORIGIN_CHOICES = (
        ("", ""),
        ("germline", "Germline"),
        ("somatic", "Somatic"),
        ("other", "Origin Other"),
    )
    allele_origin = forms.ChoiceField(choices=ALLELE_ORIGIN_CHOICES, required=False)


class ClinicalSignificanceForm(forms.Form):
    other = forms.BooleanField(required=False)
    benign = forms.BooleanField(required=False)
    likely_benign = forms.BooleanField(required=False)
    vus = forms.BooleanField(required=False)
    likely_pathogenic = forms.BooleanField(required=False)
    pathogenic = forms.BooleanField(required=False)

    helper = form_helper_horizontal()


class SampleClassificationForm(forms.Form):
    max_samples = forms.IntegerField(required=False, initial=1, min_value=1)
    max_results = forms.IntegerField(required=False, initial=5, min_value=1)
    hom_ref = forms.BooleanField(required=False)
    het = forms.BooleanField(required=False, initial=True)
    hom_alt = forms.BooleanField(required=False, initial=True)
