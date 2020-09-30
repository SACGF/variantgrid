from dal import autocomplete
from django import forms
from django.forms.formsets import BaseFormSet
from django.forms.models import ALL_FIELDS, ModelForm
from django.forms.widgets import TextInput

from library.forms import ROFormMixin
from pedigree.models import PedFile, PedFileFamily, Pedigree
from snpdb.models import Cohort


class ReadOnlyPedFileForm(ModelForm, ROFormMixin):
    class Meta:
        model = PedFile
        fields = ALL_FIELDS
        read_only = ('user', 'import_status')
        widgets = {'name': TextInput()}


class PedigreeForm(ModelForm, ROFormMixin):
    class Meta:
        model = Pedigree
        fields = ALL_FIELDS
        read_only = ('user', )
        widgets = {'name': TextInput()}


class UserCohortandPedFileFamilyForm(forms.Form):
    ped_file_family = forms.ModelChoiceField(queryset=PedFileFamily.objects.all())
    cohort = forms.ModelChoiceField(queryset=Cohort.objects.all(),
                                    widget=autocomplete.ModelSelect2(url='cohort_autocomplete',
                                                                     attrs={'data-placeholder': 'Cohort...'}))

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        self.fields['cohort'].queryset = Cohort.filter_for_user(user)
        self.fields['ped_file_family'].queryset = PedFileFamily.filter_for_user(user)


class BaseCohortSamplesForPedFileRecordsFormSet(BaseFormSet):

    def clean(self):
        """ Checks that no two records have the same cohort sample """
        if any(self.errors):
            # Don't bother validating the formset unless each form is valid on its own
            return
        cohort_samples = set()
        for form in self.forms:
            cs = form.cleaned_data['cohort_sample']
            if cs is not None:
                if cs in cohort_samples:
                    raise forms.ValidationError(f"Can only assign each cohort sample once (duplicate: {cs})")
                cohort_samples.add(cs)
