from dal import forward
from django import forms
from django.forms.models import inlineformset_factory, modelform_factory, modelformset_factory
from django.forms.widgets import TextInput

from library.django_utils.autocomplete_utils import ModelSelect2
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import (
    Extraction,
    ExternalPK,
    Patient,
    PatientModification,
    PatientRecordOriginType,
    Specimen,
)
from patients.models_enums import PopulationGroup


class PatientForm(forms.ModelForm):
    population = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
        choices=PopulationGroup.choices,
    )

    class Meta:
        model = Patient
        # patient_code first - the de-identified code we want people to enter and that's shown everywhere
        fields = ['patient_code', 'family_code', 'first_name', 'last_name',
                  'date_of_birth', 'date_of_death', 'sex',
                  'consanguineous', 'affected', 'phenotype']
        widgets = {'first_name': TextInput(attrs={'placeholder': 'First Name'}),
                   'last_name': TextInput(attrs={'placeholder': 'Last Name'}),
                   'family_code': TextInput(attrs={'placeholder': 'Family Code'}),
                   'patient_code': TextInput(attrs={'placeholder': 'Patient Code (de-identified)'}),
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
            created = patient.pk is None
            patient.save(phenotype_approval_user=self.user)
            if created:
                assign_permission_to_user_and_groups(self.user, patient)

            pop_set = patient.patientpopulation_set
            old_populations = set(pop_set.all().values_list("population", flat=True))

            changed = []
            for f in self.changed_data:
                if f in self.old_patient_data:  # Only care about fields on the patient model
                    old_val = self.old_patient_data[f]
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


# Tissue has no way to be created yet (#1747), so both editors here leave out a dropdown that
# could only ever be empty
PatientSpecimenFormSet = inlineformset_factory(Patient,
                                               Specimen,
                                               can_delete=True,
                                               exclude=['external_pk', 'tissue'],
                                               widgets={'name': TextInput(),
                                                        'description': TextInput(),
                                                        'reference_id': TextInput(),
                                                        'collected_by': TextInput(),
                                                        'collection_date': TextInput(attrs={'class': 'date-picker'}),
                                                        'received_date': TextInput(attrs={'class': 'date-picker'})},
                                               extra=1)


# The specimen/extraction pages edit the same fields as the patient tabs' formsets - the tabs stay the
# bulk editor, external_pk is the tracking system's identity so both show it read only beside the form
SpecimenForm = modelform_factory(Specimen,
                                 exclude=['external_pk', 'patient', 'tissue'],
                                 widgets={'description': TextInput(),
                                          'reference_id': TextInput(),
                                          'collected_by': TextInput(),
                                          'collection_date': TextInput(attrs={'class': 'date-picker'}),
                                          'received_date': TextInput(attrs={'class': 'date-picker'})})


class ExtractionForm(forms.ModelForm):
    """ The specimen page and the extraction page both work on a specimen you already have in hand, so
        it isn't a field here - which is why the unique_together it takes part in is checked below """

    class Meta:
        model = Extraction
        fields = ['reference_id', 'nucleic_acid_source', 'extraction_date']
        widgets = {'reference_id': TextInput(),
                   'extraction_date': TextInput(attrs={'class': 'date-picker'})}

    def clean_reference_id(self):
        # Unnamed extractions are all distinct under Postgres (@see Extraction.Meta), so an empty box
        # means unnamed rather than named ""
        return self.cleaned_data['reference_id'] or None

    def clean(self):
        cleaned_data = super().clean()
        # Django skips a unique_together whose fields the form leaves out, so a clash would otherwise
        # only surface as an IntegrityError
        reference_id = cleaned_data.get('reference_id')
        if reference_id and self.instance.specimen_id:
            clash = Extraction.objects.filter(specimen_id=self.instance.specimen_id,
                                              reference_id=reference_id).exclude(pk=self.instance.pk)
            if clash.exists():
                self.add_error('reference_id', "This specimen already has an extraction with this reference.")
        return cleaned_data


def patient_extraction_formset_factory(patient):
    """ Not an inline formset as Extraction hangs off Specimen rather than Patient, so the patient is
        forwarded as a constant to narrow the specimen autocomplete. The view sets the queryset,
        which is what actually restricts what can be saved """
    specimen_widget = ModelSelect2(url='specimen_autocomplete',
                                   forward=(forward.Const(patient.pk, 'patient'),),
                                   attrs={'data-placeholder': 'Specimen...'})
    return modelformset_factory(Extraction,
                                form=ExtractionForm,
                                can_delete=True,
                                fields=['specimen', 'reference_id', 'nucleic_acid_source', 'extraction_date'],
                                widgets={'specimen': specimen_widget,
                                         'reference_id': TextInput(),
                                         'extraction_date': TextInput(attrs={'class': 'date-picker'})},
                                extra=1)


def external_pk_autocomplete_form_factory(external_type):
    class ExternalPKNameForm(forms.Form):
        external_pk = forms.ModelChoiceField(queryset=ExternalPK.objects.all(),
                                             widget=ModelSelect2(url='external_pk_autocomplete',
                                                                 attrs={'data-placeholder': f"{external_type}..."},
                                                                 forward=(forward.Const(external_type, 'external_type'),)))

    return ExternalPKNameForm()
