from django import forms
from django.forms.widgets import HiddenInput

from genes.models import GeneSymbol
from library.django_utils.autocomplete_utils import ModelSelect2Multiple


class AllVariantsGeneSymbolForm(forms.Form):
    """ The All Variants page's gene filter - the widget provides the autocomplete, the page reads the
        selection straight out of the DOM """
    gene_symbols = forms.ModelMultipleChoiceField(
        required=False,
        label="Gene Symbols",
        queryset=GeneSymbol.objects.all(),
        widget=ModelSelect2Multiple(url='gene_symbol_autocomplete',
                                    attrs={'data-placeholder': 'Gene Symbol...'}))


class SearchForm(forms.Form):
    search = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'Search...'}), required=True)
    mode = forms.CharField(required=False)

    def __init__(self, *args, **kwargs):
        search_allow_blank = kwargs.pop("search_allow_blank", False)
        super().__init__(*args, **kwargs)
        if search_allow_blank:
            self.fields['search'].required = False

    @property
    def classify(self):
        return False


class SearchAndClassifyForm(forms.Form):
    search = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'Search...'}), required=True)
    classify = forms.BooleanField(widget=HiddenInput(), required=False)

    def __init__(self, *args, **kwargs):
        search_allow_blank = kwargs.pop("search_allow_blank", False)
        super().__init__(*args, **kwargs)
        if search_allow_blank:
            self.fields['search'].required = False
