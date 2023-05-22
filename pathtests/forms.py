from dal import forward
from django import forms
from django.core.exceptions import ValidationError
from django.forms.widgets import TextInput

from genes.models import GeneList
from library.django_utils.autocomplete_utils import ModelSelect2
from library.utils import is_not_none
from pathtests.models import PathologyTest, PathologyTestVersion, Case, \
    PathologyTestOrder
from seqauto.models import EnrichmentKit


class ActivePathologyTestForm(forms.Form):
    pathology_test = forms.ModelChoiceField(queryset=PathologyTest.objects.all(),
                                            widget=ModelSelect2(url='pathology_test_autocomplete',
                                                                attrs={'data-placeholder': 'Pathology Test...'},
                                                                forward=(forward.Const(True, 'active'),)))


class SelectPathologyTestForm(forms.Form):
    pathology_test = forms.ModelChoiceField(queryset=PathologyTest.objects.all(),
                                            widget=ModelSelect2(url='pathology_test_autocomplete',
                                                                attrs={'data-placeholder': 'Pathology Test...'}))


class SelectPathologyTestVersionForm(forms.Form):
    pathology_test_version = forms.ModelChoiceField(queryset=PathologyTestVersion.objects.all(),
                                                    widget=ModelSelect2(url='pathology_test_version_autocomplete',
                                                                        attrs={'data-placeholder': 'Pathology Test Version...'}))


class PathologyTestVersionForm(forms.ModelForm):
    enrichment_kit = forms.ModelChoiceField(queryset=EnrichmentKit.objects.all(),
                                            widget=ModelSelect2(url='enrichment_kit_autocomplete',
                                                                attrs={'data-placeholder': 'Enrichment Kit...'}))

    class Meta:
        model = PathologyTestVersion
        fields = ('enrichment_kit',)


class CreatePathologyTestForm(forms.Form):
    name = forms.CharField(required=True)
    gene_list = forms.ModelChoiceField(queryset=GeneList.objects.all(),
                                       required=False,
                                       widget=ModelSelect2(url='gene_list_autocomplete',
                                                           attrs={'data-placeholder': 'Gene List...'}), )
    pathology_test_version = forms.ModelChoiceField(queryset=PathologyTestVersion.objects.all(),
                                                    required=False,
                                                    widget=ModelSelect2(url='pathology_test_version_autocomplete',
                                                                        attrs={
                                                                            'data-placeholder': 'Pathology Test Version...'}))

    def clean_name(self):
        name = self.cleaned_data['name']
        if PathologyTest.objects.filter(name=name).exists():
            msg = f"An existing Pathology Test called '{name}' already exists"
            raise ValidationError(msg)
        return name

    def clean(self):
        cleaned_data = super().clean()

        EXCLUSIVE_FIELDS = ['gene_list', 'pathology_test_version']
        field_values = [self.cleaned_data[s] for s in EXCLUSIVE_FIELDS]
        if len(list(filter(is_not_none, field_values))) > 1:
            exclusive_fields = ', '.join(EXCLUSIVE_FIELDS)
            raise ValidationError(f"You must select at most ONE of {exclusive_fields}")
        return cleaned_data


class PathologyTestOrderForm(forms.ModelForm):
    class Meta:
        model = PathologyTestOrder
        fields = '__all__'
        widgets = {'name': TextInput(),
                   'external_pk': ModelSelect2(url='external_pk_autocomplete',
                                               attrs={'data-placeholder': 'External ID...'}),
                   'custom_gene_list': ModelSelect2(url='gene_autocomplete',
                                                    attrs={'data-placeholder': 'Gene List...'}),
                   'pathology_test_version': ModelSelect2(url='pathology_test_version_autocomplete',
                                                          attrs={'data-placeholder': 'Pathology Test Version...'}),
                   'user': ModelSelect2(url='user_autocomplete',
                                        attrs={'data-placeholder': 'User...'})}


class CaseForm(forms.ModelForm):
    class Meta:
        model = Case
        fields = '__all__'
        widgets = {'name': TextInput(),
                   'lead_scientist': ModelSelect2(url='user_autocomplete',
                                                  attrs={'data-placeholder': 'User...'}),
                   'patient': ModelSelect2(url='user_autocomplete',
                                           attrs={'data-placeholder': 'Patient...'})}

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        if not self.instance.can_write(user):
            for f in self.fields.values():
                f.disabled = True
