from dal import forward
from django import forms
from django.forms.models import inlineformset_factory, ALL_FIELDS
from django.forms.widgets import TextInput

from library.django_utils.autocomplete_utils import ModelSelect2
from patients.models import Patient, Specimen, ExternalPK, PatientModification, \
    PatientRecordOriginType
from patients.models_enums import PopulationGroup


class PatientForm(forms.ModelForm):
    population = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
        choices=PopulationGroup.choices,
    )

    class Meta:
        model = Patient
        fields = ['first_name', 'last_name', 'family_code',
                  'date_of_birth', 'date_of_death', 'sex',
                  'consanguineous', 'affected', 'phenotype']
        widgets = {'first_name': TextInput(attrs={'placeholder': 'First Name'}),
                   'last_name': TextInput(attrs={'placeholder': 'Last Name'}),
                   'family_code': TextInput(attrs={'placeholder': 'Family Code'}),
                   'date_of_birth': TextInput(attrs={'class': 'date-picker', 'placeholder': 'Date of Birth'}),
                   'date_of_death': TextInput(attrs={'class': 'date-picker', 'placeholder': 'Date of Death'})}

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        instance = kwargs.get("instance")
        initial = kwargs.get("initial") or {}
        kwargs["initial"] = initial
        kwargs.setdefault('label_suffix', '')
        if instance and instance.pk:
            self.old_patient_data = instance.__dict__.copy()
            pop_set = instance.patientpopulation_set.all()
            initial["population"] = list(set(pop_set.values_list("population", flat=True)))
        else:
            self.old_patient_data = {}

        super().__init__(*args, **kwargs)

    def create_patient_modification(self, description):
        PatientModification.objects.create(patient=self.instance,
                                           user=self.user,
                                           description=description,
                                           origin=PatientRecordOriginType.MANUAL_VG_GUI)

    def save(self, commit=True):
        patient = super().save(commit=False)
        if commit:
            patient.save(phenotype_approval_user=self.user)
            pop_set = patient.patientpopulation_set
            old_populations = set(pop_set.all().values_list("population", flat=True))

            changed = []
            for f in self.changed_data:
                old_val = self.old_patient_data.get(f)
                if old_val:  # Only care about fields on the patient model
                    new_val = getattr(patient, f)
                    changed.append(f"{f}: '{old_val}' to '{new_val}'")

            if changed:
                description = "Changed: " + ', '.join(changed)
                self.create_patient_modification(description)

            new_populations = set(self.cleaned_data['population'])
            if new_populations != old_populations:
                added = new_populations - old_populations
                removed = old_populations - new_populations

                population_changes = "Population - "
                if added:
                    population_changes += "added: " + ",".join(added)
                    for population in added:
                        pop_set.create(population=population)

                if removed:
                    population_changes += " removed: " + ",".join(removed)
                    for population in removed:
                        pop_set.filter(population=population).delete()

                if self.old_patient_data:  # Don't save changes if newly created
                    self.create_patient_modification(population_changes)

        return patient


class PatientContactForm(forms.ModelForm):
    class Meta:
        model = Patient
        fields = ['street_address', 'suburb', 'postcode', 'telephone']
        widgets = {'street_address': TextInput(attrs={'placeholder': 'Number, Street'}),
                   'suburb': TextInput(attrs={'placeholder': 'Suburb'}),
                   'postcode': TextInput(attrs={'placeholder': 'Postcode'}),
                   'telephone': TextInput(attrs={'placeholder': 'Numbers only'})}


class PatientSearchForm(forms.Form):
    patient = forms.ModelChoiceField(queryset=Patient.objects.all(),
                                     widget=ModelSelect2(url='patient_autocomplete',
                                                         attrs={'data-placeholder': 'Patient...'}))
    family_code = forms.CharField(widget=TextInput(attrs={'placeholder': 'Family Code'}))
    phenotype = forms.CharField(widget=TextInput(attrs={'placeholder': 'Phenotype text'}))


PatientSpecimenFormSet = inlineformset_factory(Patient,
                                               Specimen,
                                               can_delete=True,
                                               fields=ALL_FIELDS,
                                               widgets={'name': TextInput(),
                                                        'description': TextInput(),
                                                        'reference_id': TextInput(),
                                                        'collected_by': TextInput(),
                                                        'collection_date': TextInput(attrs={'class': 'date-picker'}),
                                                        'received_date': TextInput(attrs={'class': 'date-picker'})},
                                               extra=1)


def external_pk_autocomplete_form_factory(external_type):
    class ExternalPKNameForm(forms.Form):
        external_pk = forms.ModelChoiceField(queryset=ExternalPK.objects.all(),
                                             widget=ModelSelect2(url='external_pk_autocomplete',
                                                                 attrs={'data-placeholder': f"{external_type}..."},
                                                                 forward=(forward.Const(external_type, 'external_type'),)))

    return ExternalPKNameForm()
