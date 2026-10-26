from django import forms
from django.db.models import Q
from django.forms.formsets import BaseFormSet
from django.forms.models import ALL_FIELDS, ModelForm
from django.forms.widgets import TextInput

from library.django_utils.autocomplete_utils import ModelSelect2
from library.forms import ROFormMixin
from pedigree.models import PedFile, PedFileFamily, Pedigree
from snpdb.models import Cohort


class ReadOnlyPedFileForm(ModelForm, ROFormMixin):
    class Meta:
        model = PedFile
        fields = ['name', 'user']
        read_only_display = ('user', )
        widgets = {'name': TextInput()}


class PedigreeForm(ModelForm, ROFormMixin):
    class Meta:
        model = Pedigree
        fields = ALL_FIELDS
        read_only_display = ('user', )
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        cohort_qs = Cohort.filter_for_user(user)
        ped_file_family_qs = PedFileFamily.filter_for_user(user)
        if self.instance.pk:
            # Sharing a pedigree doesn't share its cohort or PED file, so keep the current ones selectable
            cohort_qs = Cohort.objects.filter(Q(pk__in=cohort_qs.values("pk")) | Q(pk=self.instance.cohort_id))
            ped_file_family_qs = PedFileFamily.objects.filter(Q(pk__in=ped_file_family_qs.values("pk"))
                                                              | Q(pk=self.instance.ped_file_family_id))
        self.fields['cohort'].queryset = cohort_qs
        self.fields['ped_file_family'].queryset = ped_file_family_qs


class UserCohortandPedFileFamilyForm(forms.Form):
    ped_file_family = forms.ModelChoiceField(queryset=PedFileFamily.objects.all())
    cohort = forms.ModelChoiceField(queryset=Cohort.objects.all(),
                                    widget=ModelSelect2(url='cohort_autocomplete',
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
