from crispy_forms.helper import FormHelper
from django import forms

from library.django_utils.autocomplete_utils import ModelSelect2, ModelSelect2Multiple
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


class MONDOForm(forms.Form):
    mondo = forms.ModelChoiceField(queryset=OntologyTerm.objects.all(),
                                   required=False,
                                   widget=ModelSelect2(url='mondo_autocomplete',
                                                       attrs={'data-placeholder': 'MONDO...'}))


class HGNCForm(forms.Form):
    hgnc = forms.ModelChoiceField(queryset=OntologyTerm.objects.all(),
                                  required=False,
                                  widget=ModelSelect2(url='hgnc_autocomplete',
                                                      attrs={'data-placeholder': 'Gene/HGNC...'}))


class PhenotypeMultipleSelectForm(forms.Form):
    omim = forms.ModelMultipleChoiceField(required=False,
                                          queryset=OntologyTerm.objects.all(),
                                          widget=ModelSelect2Multiple(url='omim_autocomplete',
                                                                      attrs={'data-placeholder': 'OMIM...',
                                                                             'class': 'omim'}))
    hpo = forms.ModelMultipleChoiceField(required=False,
                                         queryset=OntologyTerm.objects.all(),
                                         widget=ModelSelect2Multiple(url='hpo_autocomplete',
                                                                     attrs={'data-placeholder': 'HPO...',
                                                                            'class': 'hpo'}))

    mondo = forms.ModelMultipleChoiceField(required=False,
                                           queryset=OntologyTerm.objects.all(),
                                           widget=ModelSelect2Multiple(url='mondo_autocomplete',
                                                                       attrs={'data-placeholder': 'MONDO...',
                                                                              'class': 'mondo'}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_class = 'form-horizontal'
        self.helper.label_class = 'col-sm-2'  # Bootstrap column for labels
        self.helper.field_class = 'col-sm-10'  # Bootstrap column for inputs
