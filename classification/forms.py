from django import forms
from django.forms.widgets import TextInput

from classification.models.classification_settings_models import LabClassificationSettings
from library.forms import ROFormMixin
from uicore.utils.form_helpers import form_helper_horizontal


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
    pop_pm2_max_af = forms.FloatField(required=False, initial=3e-4)

    # ClinVar filters
    clinvar_min_conflict_distance = forms.IntegerField(required=False, initial=2)
    clinvar_min_stars = forms.IntegerField(required=False, initial=2)

    computational_vus_spliceai_min = forms.FloatField(required=False, initial=0.5)


class LabClassificationSettingsForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = LabClassificationSettings
        fields = ("consolidates_variant_classifications", "upload_location", "upload_automatic", "upload_instructions")
        read_only = ("upload_location", "upload_automatic", "upload_instructions")
        widgets = {
            "upload_location": TextInput(),
        }
        labels = {
            "consolidates_variant_classifications": "Consolidates Variant Classifications",
            "upload_location": "Upload Location",
            "upload_automatic": "Upload Automatic",
            "upload_instructions": "Upload Instructions",
        }
        help_texts = {
            "consolidates_variant_classifications": "The lab makes one classification per variant rather than one per patient, so patient counts are a minimum.",
            "upload_location": "If provided, classification uploads can be done via the classifications/upload page.",
            "upload_automatic": "Uploaded files are processed automatically, otherwise they wait for manual review.",
        }
