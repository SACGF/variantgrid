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


class ClassificationEvidenceUpdateForm(forms.Form):
    max_results = forms.IntegerField(required=False, initial=10, min_value=1)
    population = forms.BooleanField(required=False, initial=True)
    clinvar = forms.BooleanField(required=False, initial=True)
    computational = forms.BooleanField(required=False, initial=True)
    # gene_disease = forms.BooleanField(required=False, initial=True)

    # Pop
    pop_no_ba1_min_af = forms.FloatField(required=False, initial=0.01)
    pop_no_bs1_min_af = forms.FloatField(required=False, initial=0.001)
    pop_recessive_no_bs2_min_homozygotes = forms.IntegerField(required=False, initial=1)
    pop_pm2_min_af = forms.FloatField(required=False, initial=1e-4)

    # ClinVar filters
    clinvar_min_conflict_distance = forms.IntegerField(required=False, initial=2)
    clinvar_min_stars = forms.IntegerField(required=False, initial=2)

    computational_vus_spliceai_min = forms.FloatField(required=False, initial=0.5)

