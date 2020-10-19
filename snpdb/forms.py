import collections

from crispy_forms.bootstrap import FieldWithButtons
from crispy_forms.layout import Layout, Submit, Field
from dal import autocomplete, forward
from django import forms
from django.conf import settings
from django.contrib.auth.models import User
from django.forms.forms import DeclarativeFieldsMetaclass
from django.forms.widgets import TextInput, HiddenInput, NullBooleanSelect
from guardian import shortcuts
from guardian.shortcuts import assign_perm, remove_perm

from library.forms import ROFormMixin
from library.guardian_utils import DjangoPermission
from snpdb import models
from snpdb.models import VCF, Sample, Cohort, VARIANT_PATTERN, \
    DBSNP_PATTERN, UserContact, Tag, UserSettings, GenomicIntervalsCollection, ImportStatus, \
    SettingsInitialGroupPermission, LabUserSettingsOverride, UserSettingsOverride, \
    OrganizationUserSettingsOverride
from snpdb.models.models import Lab, Organization
from snpdb.models.models_genome import GenomeBuild
from uicore.utils.form_helpers import form_helper_horizontal, FormHelperHelper
from variantgrid.perm_path import get_visible_url_names


class BaseForm(forms.BaseForm):
    """
    Removes the ":" suffix on forms
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('label_suffix', '')
        super().__init__(*args, **kwargs)


class BaseDeclareForm(BaseForm, metaclass=DeclarativeFieldsMetaclass):
    pass


class BaseModelForm(forms.ModelForm):
    """
    Removes the ":" suffix on forms
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('label_suffix', '')
        super().__init__(*args, **kwargs)


class GenomeBuildAutocompleteForwardMixin:
    genome_build_fields = []

    def __init__(self, *args, **kwargs):
        genome_build = kwargs.pop("genome_build", None)
        super().__init__(*args, **kwargs)
        if genome_build:
            for f in self.genome_build_fields:
                self.fields[f].widget.forward.append(forward.Const(genome_build.pk, "genome_build_id"))


class UserSelectForm(forms.Form):
    user = forms.ModelChoiceField(queryset=User.objects.all(),
                                  required=False,
                                  widget=autocomplete.ModelSelect2(url='username_autocomplete',
                                                                   attrs={'data-placeholder': 'Username...'}))


class KeycloakUserForm(BaseForm):
    base_fields = collections.OrderedDict()
    base_fields['first_name'] = forms.CharField(required=True, max_length=100)
    base_fields['last_name'] = forms.CharField(required=True, max_length=100)
    base_fields['email'] = forms.EmailField(required=True)
    base_fields['lab'] = forms.ModelChoiceField(
        required=True,
        queryset=Lab.objects.all(),
        widget=autocomplete.ModelSelect2(url='lab_autocomplete',
                                         attrs={'data-placeholder': 'Lab...'}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = form_helper_horizontal()
        self.helper.layout = Layout(
            Field('first_name'),
            Field('last_name'),
            Field('email'),
            Field('lab'),
            Submit('submit', 'Create User')
        )


class LabSelectForm(forms.Form):
    lab = forms.ModelChoiceField(queryset=Lab.objects.all(),
                                 required=False,
                                 widget=autocomplete.ModelSelect2(url='lab_autocomplete',
                                                                  attrs={'data-placeholder': 'Lab...'}))


class TagForm(forms.Form):
    tag = forms.ModelChoiceField(queryset=Tag.objects.all(),
                                 widget=autocomplete.ModelSelect2(url='tag_autocomplete',
                                                                  attrs={'data-placeholder': 'Tag...'}))


class UserForm(BaseModelForm):

    class Meta:
        model = User
        fields = ['first_name', 'last_name', 'email']


class UserContactForm(BaseModelForm):

    class Meta:
        model = UserContact
        fields = ['phone_number']
        widgets = {'phone_number': TextInput()}


class LabForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = Lab
        # fields = '__all__'
        exclude = ("classification_config", )
        read_only = ("name", "external", "group_name", "organization", "upload_location")
        widgets = {
            "name": TextInput(),
            "city": TextInput(),
            "state": TextInput(),
            "country": TextInput(),
            "url": TextInput(),
            "css_class": TextInput(),
            "group_name": TextInput(),
            "upload_location": TextInput()
        }


class OrganizationForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = Organization
        fields = '__all__'
        read_only = ("name", "group_name", "active")
        widgets = {
            "name": TextInput(),
            "short_name": TextInput(),
            "group_name": TextInput(),
        }


class VCFForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = models.VCF
        fields = ['name', 'date', 'genome_build', 'user', 'project', 'import_status']
        read_only = ('import_status', )
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.instance.genome_build:
            self.fields["genome_build"].disabled = True


class GroupPermissionForm(forms.Form):
    """ Displays and allows changing Guardian permissions for an object """
    read = forms.BooleanField(required=False)
    write = forms.BooleanField(required=False)

    def __init__(self, *args, **kwargs):
        self.group = kwargs.pop("group")
        obj = kwargs.pop("obj")

        self.obj = obj
        self.read_perm_name = DjangoPermission.perm(obj, DjangoPermission.READ)
        self.write_perm_name = DjangoPermission.perm(obj, DjangoPermission.WRITE)

        if self.group:
            kwargs['prefix'] = f"{self.group.name}_{self.obj.pk}"

        super().__init__(*args, **kwargs)
        self.set_checked()

    def set_checked(self):
        perms = shortcuts.get_perms(self.group, self.obj)
        read_checked = self.read_perm_name in perms
        write_checked = self.write_perm_name in perms

        self.set_field_checked('read', read_checked)
        self.set_field_checked('write', write_checked)

    def save(self):
        read = self.cleaned_data['read']
        write = self.cleaned_data['write']

        self.set_perm(self.read_perm_name, read, self.group, self.obj)
        self.set_perm(self.write_perm_name, write, self.group, self.obj)
        self.set_checked()

    @staticmethod
    def set_perm(perm, test, group_or_user, obj):
        if test:
            assign_perm(perm, group_or_user, obj)
        else:
            remove_perm(perm, group_or_user, obj)

    def set_field_checked(self, field, test):
        CHECKED = 'checked'
        attrs = self.fields[field].widget.attrs
        if test:
            attrs[CHECKED] = 'checked'
        elif attrs.get(CHECKED):
            del attrs[CHECKED]


class SettingsInitialGroupPermissionForm(forms.Form):
    """ Not done as a modelform as we always want to display all groups,
        and not make objects unless required """
    read = forms.BooleanField(required=False)
    write = forms.BooleanField(required=False)

    def __init__(self, *args, **kwargs):
        self.settings_override = kwargs.pop("settings_override")
        self.group = kwargs.pop("group")

        initial = kwargs.get("initial", {})
        self._original_read = initial.get("read")
        self._original_write = initial.get("write")

        kwargs['prefix'] = f"{self.settings_override.pk}_{self.group.name}"
        super().__init__(*args, **kwargs)

    def save(self):
        read = self.cleaned_data['read']
        write = self.cleaned_data['write']

        if read != self._original_read or write != self._original_write:
            defaults = {"read": read, "write": write}
            SettingsInitialGroupPermission.objects.update_or_create(settings=self.settings_override,
                                                                    group=self.group, defaults=defaults)


class SampleForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = models.Sample
        exclude = ['vcf', 'has_genotype']
        read_only = ('vcf_sample_name', 'import_status')
        widgets = {'vcf_sample_name': TextInput(),
                   'name': TextInput(),
                   'bam_file_path': TextInput(),
                   'patient': autocomplete.ModelSelect2(url='patient_autocomplete',
                                                        attrs={'data-placeholder': 'Patient...'}),
                   'specimen': autocomplete.ModelSelect2(url='specimen_autocomplete',
                                                         forward=['patient'],
                                                         attrs={'data-placeholder': 'Specimen...'})}

    def clean(self):
        cleaned_data = super().clean()
        specimen = cleaned_data.get('specimen')
        if specimen:
            patient = cleaned_data.get('patient')
            if patient is None:
                self.add_error('patient', "Patient must be supplied if specimen supplied")
            elif specimen.patient != patient:
                msg = "Specimen must be from supplied patient"
                self.add_error('specimen', msg)
                self.add_error('patient', msg)
        return cleaned_data


class SampleChoiceForm(GenomeBuildAutocompleteForwardMixin, BaseDeclareForm):
    genome_build_fields = ["sample"]
    sample = forms.ModelChoiceField(queryset=Sample.objects.all(),
                                    widget=autocomplete.ModelSelect2(url='sample_autocomplete',
                                                                     attrs={'data-placeholder': 'Sample...'}))


class VCFChoiceForm(GenomeBuildAutocompleteForwardMixin, BaseDeclareForm):
    genome_build_fields = ["vcf"]

    vcf = forms.ModelChoiceField(queryset=VCF.objects.all(),
                                 widget=autocomplete.ModelSelect2(url='vcf_autocomplete',
                                                                  attrs={'data-placeholder': 'VCF...'}))


class BlankNullBooleanSelect(NullBooleanSelect):
    """
    Has "-" instead of "Unknown"
    """
    def __init__(self, attrs=None):
        super().__init__(attrs)
        tup = self.choices[0]
        self.choices[0] = (tup[0], "-")


class SettingsOverrideForm(BaseModelForm):
    """ Warning: This hides some fields (in _hide_unused_fields) so not all fields from model will exist """
    genome_builds = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
    )

    class Meta:
        model = models.SettingsOverride
        fields = "__all__"
        widgets = {"email_weekly_updates": BlankNullBooleanSelect(),
                   "email_discordance_updates": BlankNullBooleanSelect(),
                   "variant_link_in_analysis_opens_new_tab": BlankNullBooleanSelect(),
                   "tool_tips": BlankNullBooleanSelect(),
                   "node_sql_tab": BlankNullBooleanSelect(),
                   "import_messages": BlankNullBooleanSelect(),
                   'default_sort_by_column': autocomplete.ModelSelect2(url='custom_column_autocomplete',
                                                                       forward=['columns'],
                                                                       attrs={'data-placeholder': 'Column...'})}

    def __init__(self, *args, **kwargs):
        instance = kwargs.get("instance")
        initial = kwargs.get("initial") or {}
        kwargs["initial"] = initial
        if instance and instance.pk:
            usgb_set = instance.settingsgenomebuild_set
            genome_builds = set(usgb_set.all().values_list("genome_build", flat=True))
            initial["genome_builds"] = list(genome_builds)

        super().__init__(*args, **kwargs)
        self.fields['genome_builds'].choices = GenomeBuild.get_choices()
        self.fields['columns'].queryset = models.CustomColumnsCollection.objects.filter(user__isnull=True)  # public
        self.fields['default_genome_build'].queryset = GenomeBuild.builds_with_annotation()
        self._hide_unused_fields()

    def _hide_unused_fields(self):
        visible_url_names = get_visible_url_names()
        analysis_enabled = visible_url_names.get("analysis")
        upload_enabled = visible_url_names.get("upload") and analysis_enabled
        igv_links_enabled = any((analysis_enabled,
                                 settings.VARIANT_DETAILS_SHOW_ANNOTATION,
                                 settings.VARIANT_DETAILS_SHOW_SAMPLES))
        field_visibility = {
            "email_weekly_updates": settings.DISCORDANCE_ENABLED,
            "email_discordance_updates": settings.DISCORDANCE_ENABLED,
            "columns": analysis_enabled,
            "default_sort_by_column": analysis_enabled,
            "variant_link_in_analysis_opens_new_tab": analysis_enabled,
            "tool_tips": analysis_enabled,
            "node_sql_tab": analysis_enabled,
            "import_messages": upload_enabled,
            "igv_port": igv_links_enabled,
        }

        for f, visible in field_visibility.items():
            if f in self.fields and not visible:
                del self.fields[f]

    def save(self, commit=True):
        settings_override = super().save(commit=commit)
        if commit:
            usgb_set = settings_override.settingsgenomebuild_set
            old_genome_builds = set(usgb_set.all().values_list("genome_build", flat=True))
            new_genome_builds = set(self.cleaned_data['genome_builds'])

            if settings_override.default_genome_build:
                # Always use default if set
                new_genome_builds.add(settings_override.default_genome_build.name)

            if new_genome_builds != old_genome_builds:
                added = new_genome_builds - old_genome_builds
                removed = old_genome_builds - new_genome_builds

                if added:
                    for build_name in added:
                        genome_build = GenomeBuild.objects.get(name=build_name)
                        usgb_set.create(genome_build=genome_build)

                if removed:
                    for build_name in removed:
                        usgb_set.filter(genome_build__name=build_name).delete()

        return settings_override


class OrganizationUserSettingsOverrideForm(SettingsOverrideForm):
    class Meta(SettingsOverrideForm.Meta):
        model = OrganizationUserSettingsOverride
        exclude = ['organization']


class LabUserSettingsOverrideForm(SettingsOverrideForm):
    class Meta(SettingsOverrideForm.Meta):
        model = LabUserSettingsOverride
        exclude = ['lab']


class UserSettingsOverrideForm(SettingsOverrideForm):
    class Meta(SettingsOverrideForm.Meta):
        model = UserSettingsOverride
        exclude = ['user', 'oauth_sub']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Some model fields may have been hidden in _hide_unused_fields - so need to check if they exist
        user = self.instance.user
        if "columns" in self.fields:
            self.fields['columns'].queryset = models.CustomColumnsCollection.filter_for_user(user)
        self.fields['default_lab'].queryset = Lab.valid_labs_qs(user)


class CreateCohortForm(BaseModelForm):

    class Meta:
        model = models.Cohort
        fields = ['name', "genome_build"]
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['genome_build'].choices = GenomeBuild.get_choices()


class CohortForm(BaseModelForm):

    class Meta:
        model = models.Cohort
        fields = ['name']
        widgets = {'name': TextInput()}


class CustomColumnsCollectionForm(forms.ModelForm):

    class Meta:
        model = models.CustomColumnsCollection
        exclude = ['user', 'version_id']
        widgets = {'name': TextInput()}

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        helper = FormHelperHelper.instance.horizontal
        self.helper = helper

        helper.layout = Layout(
            FieldWithButtons('name', Submit(name="Create", value="create", css_class="btn btn-primary"))
        )

    def save(self, commit=True):
        instance = super().save(commit=False)
        instance.user = self.user

        if commit:
            instance.save()

        return instance


class ManualVariantEntryForm(forms.Form):
    variants_text = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Variants...'}), required=True)
    genome_build = forms.ChoiceField(widget=forms.RadioSelect)

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        user_settings = UserSettings.get_for_user(user)
        initial = {"genome_build": user_settings.default_genome_build.pk}
        passed_initial = kwargs.get("initial", {})
        initial.update(passed_initial)  # overwrite
        kwargs["initial"] = initial
        super().__init__(*args, **kwargs)
        choices = GenomeBuild.get_choices()
        self.fields['genome_build'].choices = choices

    def clean(self):
        cleaned_data = super().clean()
        variants_text = cleaned_data.get("variants_text")
        if variants_text:
            for line in variants_text.split('\n'):
                line = line.strip()
                if not (DBSNP_PATTERN.match(line) or VARIANT_PATTERN.match(line)):
                    msg = f"'{line}' does not match dbSNP or variant patterns"
                    if " " in line:
                        msg += " (remember to only enter one variant per line)"
                    raise forms.ValidationError(msg)

        return cleaned_data


class UserCohortForm(forms.Form):
    cohort = forms.ModelChoiceField(queryset=Cohort.objects.all(),
                                    widget=autocomplete.ModelSelect2(url='cohort_autocomplete',
                                                                     attrs={'data-placeholder': 'Cohort...'}))


class CreateTagForm(forms.Form):
    tag = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'New Tag Name...'}), required=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        helper = FormHelperHelper.instance.horizontal
        helper.layout = Layout(
            FieldWithButtons('tag', Submit(name="Create", value="create", css_class="btn btn-primary"))
        )
        self.helper = helper


class UserSettingsGenomeBuildMixin:
    """ Mixin with ModelForm to have genome_build initialise from UserSettings """
    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        user_settings = UserSettings.get_for_user(self.user)
        initial = kwargs.get("initial", {})
        initial['genome_build'] = user_settings.default_genome_build
        kwargs["initial"] = initial
        super().__init__(*args, **kwargs)

        genome_builds_qs = user_settings.get_genome_builds()
        build_field = self.fields['genome_build']
        build_field.queryset = genome_builds_qs
        if genome_builds_qs.count() < 2:
            build_field.widget = HiddenInput()


class GenomicIntervalsCollectionForm(UserSettingsGenomeBuildMixin, forms.ModelForm, ROFormMixin):
    class Meta:
        model = GenomicIntervalsCollection
        exclude = ['category']
        read_only = ('processed_file', 'processed_records', 'import_status')
        widgets = {'name': TextInput(),
                   'processed_file': TextInput()}

    def save(self, commit=True):
        instance = super().save(commit=False)

        if commit:
            # Toggle REQUIRES_USER_INPUT if genome build set or not
            if instance.import_status == ImportStatus.REQUIRES_USER_INPUT:
                if instance.genome_build is not None:
                    instance.uploadedbed.process_bed_file()
            elif instance.import_status == ImportStatus.SUCCESS:
                if instance.genome_build is None:
                    instance.import_status = ImportStatus.REQUIRES_USER_INPUT
            instance.save()

        return instance
