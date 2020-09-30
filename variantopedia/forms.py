from django import forms
from django.forms.widgets import HiddenInput


class SearchForm(forms.Form):
    search = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'Search...'}), required=True)

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
