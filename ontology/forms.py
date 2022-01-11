from django import forms

from library.django_utils.autocomplete_utils import ModelSelect2
from ontology.models import OntologyTerm


class OMIMForm(forms.Form):
    omim = forms.ModelChoiceField(queryset=OntologyTerm.objects.all(),
                                  required=False,
                                  widget=ModelSelect2(url='omim_autocomplete',
                                                      attrs={'data-placeholder': 'OMIM...'}))


class HPOForm(forms.Form):
    hpo = forms.ModelChoiceField(queryset=OntologyTerm.objects.all(),
                                 required=False,
                                 widget=ModelSelect2(url='hpo_autocomplete',
                                                     attrs={'data-placeholder': 'HPO...'}))


class HGNCForm(forms.Form):
    hgnc = forms.ModelChoiceField(queryset=OntologyTerm.objects.all(),
                                  required=False,
                                  widget=ModelSelect2(url='hgnc_autocomplete',
                                                      attrs={'data-placeholder': 'Gene/HGNC...'}))
