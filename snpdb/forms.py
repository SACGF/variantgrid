import collections
from functools import cached_property

from crispy_forms.bootstrap import FieldWithButtons
from crispy_forms.layout import Layout, Submit, Field
from dal import forward
from django import forms
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.forms import EmailInput, URLInput, inlineformset_factory, ALL_FIELDS
from django.forms.forms import DeclarativeFieldsMetaclass
from django.forms.widgets import TextInput, HiddenInput, NullBooleanSelect
from guardian import shortcuts
from guardian.shortcuts import assign_perm, remove_perm

from annotation.models import ManualVariantEntry
from annotation.models.models_enums import ManualVariantEntryType
from library.cache import timed_cache
from library.django_utils.autocomplete_utils import ModelSelect2, ModelSelect2Multiple
from library.forms import ROFormMixin
from library.guardian_utils import DjangoPermission
from snpdb import models
from snpdb.models import VCF, Sample, Cohort, UserContact, Tag, UserSettings, GenomicIntervalsCollection, \
    ImportStatus, SettingsInitialGroupPermission, LabUserSettingsOverride, UserSettingsOverride, \
    OrganizationUserSettingsOverride, CustomColumnsCollection, Project, VariantsType, SampleFilePath
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
        for f in self.genome_build_fields:
            widget_forward = []
            if genome_build:
                widget_forward.append(forward.Const(genome_build.pk, "genome_build_id"))
            self.fields[f].widget.forward = widget_forward


class UserSelectForm(forms.Form):
    user = forms.ModelChoiceField(queryset=User.objects.all(),
                                  required=False,
                                  widget=ModelSelect2(url='username_autocomplete',
                                                      attrs={'data-placeholder': 'Username...'}))


class KeycloakUserForm(BaseForm):
    base_fields = collections.OrderedDict()
    base_fields['first_name'] = forms.CharField(required=True, max_length=100)
    base_fields['last_name'] = forms.CharField(required=True, max_length=100)
    base_fields['email'] = forms.EmailField(required=True)
    base_fields['lab'] = forms.ModelChoiceField(
        required=True,
        queryset=Lab.objects.all(),
        widget=ModelSelect2(url='lab_autocomplete',
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
                                 widget=ModelSelect2(url='lab_autocomplete',
                                                     attrs={'data-placeholder': 'Lab...'}))


class LabMultiSelectForm(forms.Form):
    lab = forms.ModelMultipleChoiceField(queryset=Lab.objects.all(),
                                         required=False,
                                         widget=ModelSelect2Multiple(url='lab_autocomplete',
                                                                     attrs={'data-placeholder': 'Lab...'}))


class TagForm(forms.Form):
    tag = forms.ModelChoiceField(queryset=Tag.objects.all(),
                                 widget=ModelSelect2(url='tag_autocomplete',
                                                     attrs={'data-placeholder': 'Tag...'}))


class UserForm(BaseModelForm):

    class Meta:
        model = User
        fields = ['first_name', 'last_name', 'email']
        labels = {
            "first_name": "First Name",
            "last_name": "Last Name",
            "email": "Email Address",
            "username": "Username"
        }


class UserContactForm(BaseModelForm):

    class Meta:
        model = UserContact
        fields = ['phone_number']
        widgets = {'phone_number': TextInput()}
        labels = {
            "phone_number": "Phone Number"
        }


class LabForm(forms.ModelForm, ROFormMixin):
    class Meta:
        model = Lab
        # fields = '__all__'
        exclude = ("classification_config", "css_class")
        read_only = ("name", "external", "group_name", "organization", "upload_location", "upload_automatic", "upload_instructions", "clinvar_key")
        widgets = {
            "name": TextInput(),
            "city": TextInput(),
            "state": TextInput(),
            "country": TextInput(),
            "contact_name": TextInput(),
            "contact_email": TextInput(),
            "contact_phone": TextInput(),
            "url": TextInput(),
            "css_class": TextInput(),
            "group_name": TextInput(),
            "upload_location": TextInput(),
            "email": EmailInput(),
            "slack_webhook": URLInput()
        }
        labels = {
            "email": "Group Email",
            "url": "URL",
            "clinvar_key": "ClinVar Key",
            "group_name": "Group Name",
            "organization": "Organisation",  # really need to do translations, even if just en-US vs en-UK
            "upload_location": "Upload Location",
            "upload_auto_pattern": "Upload Auto Pattern",
            "slack_webhook": "Slack Webhook",
            "contact_name": "Contact Name",
            "contact_email": "Contact Email",
            "contact_phone": "Contact Phone"
        }
        help_texts = {
            "email": "Lab wide email for discordance and general communications.",
            "contact_name": "Name of contact person available for other labs to contact (if applicable)",
            "contact_email": "Email address to be provided to other labs when needing to communicate",
            "contact_phone": "Phone number to be provided to other labs when needing to communicate",
            "upload_location": "If provided, classification uploads can be done via the classifications/upload page.",
            "upload_auto_pattern": "If provided, then uploading files that match this pattern will be automatically processed, otherwise there will be a delay for manual review.",
            "slack_webhook": "If provided, discordance and general communications can be posted to your Slack instance. Should look like https://hooks.slack.com/services/ABC/DEF/GHI",
            "clinvar_key": "Required to submit to ClinVar. Ask the admins if your lab is ready to submit."
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
        read_only = ('date', 'import_status')
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
    genome_build = forms.CharField()  # From VCF.genome_build
    grid_sample_label = forms.CharField(help_text="Calculated from your current user settings. May be different in analysis due to analysis settings")

    class Meta:
        model = models.Sample
        exclude = ['vcf', 'has_genotype']
        read_only = ('genome_build', 'vcf_sample_name', 'import_status', 'grid_sample_label')
        widgets = {'vcf_sample_name': TextInput(),
                   'name': TextInput(),
                   'patient': ModelSelect2(url='patient_autocomplete',
                                           attrs={'data-placeholder': 'Patient...'}),
                   'specimen': ModelSelect2(url='specimen_autocomplete',
                                            forward=['patient'],
                                            attrs={'data-placeholder': 'Specimen...'})}

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        super().__init__(*args, **kwargs)
        fields = ["genome_build"] + [f for f in self.fields if f != 'genome_build']
        self.order_fields(fields)
        self.fields['genome_build'].initial = self.instance.vcf.genome_build
        self.fields['grid_sample_label'].initial = self._get_sample_label(user, self.instance)

    @staticmethod
    def _get_sample_label(user, sample):
        user_settings = UserSettings.get_for_user(user)
        grid_sample_label_template = user_settings.grid_sample_label_template
        sample_formatter = Sample._get_sample_formatter_func(grid_sample_label_template)
        return sample_formatter(sample)

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


SampleFilesFormSet = inlineformset_factory(Sample,
                                           SampleFilePath,
                                           can_delete=True,
                                           fields=ALL_FIELDS,
                                           widgets={'label': TextInput(),
                                                    'file_path': TextInput()},
                                           extra=1)


class ProjectChoiceForm(forms.Form):
    project = forms.ModelChoiceField(queryset=Project.objects.all(),
                                     widget=ModelSelect2(url='project_autocomplete',
                                                         attrs={'data-placeholder': 'Project...'}))


class SampleChoiceForm(GenomeBuildAutocompleteForwardMixin, BaseDeclareForm):
    genome_build_fields = ["sample"]
    sample = forms.ModelChoiceField(queryset=Sample.objects.all(),
                                    widget=ModelSelect2(url='sample_autocomplete',
                                                        attrs={'data-placeholder': 'Sample...'}))


class VariantsTypeMultipleChoiceForm(forms.Form):
    variants_type = forms.MultipleChoiceField(choices=VariantsType.choices, widget=forms.CheckboxSelectMultiple())


class VCFChoiceForm(GenomeBuildAutocompleteForwardMixin, BaseDeclareForm):
    genome_build_fields = ["vcf"]

    vcf = forms.ModelChoiceField(queryset=VCF.objects.all(),
                                 widget=ModelSelect2(url='vcf_autocomplete',
                                                     attrs={'data-placeholder': 'VCF...'}))


class BlankNullBooleanSelect(NullBooleanSelect):
    """
    Has "-" instead of "Unknown"
    """
    def __init__(self, attrs=None):
        super().__init__(attrs)
        tup = self.choices[0]
        self.choices[0] = (tup[0], "-")


class SettingsFormFeatures:
    """
    Calculates which features (as relevant to the settings override forms) are enabled.
    Probably should be its own independent named class as the fact that the features are referenced
    on the settings form are second to if they're enabled features or not.
    Alternatively could just get an instance of SettingsOverrideForm to see if igv_links_enabled but that seems like overkill
    """

    @cached_property
    def analysis_enabled(self) -> bool:
        return get_visible_url_names().get("analysis")

    @cached_property
    def upload_enabled(self) -> bool:
        return self.analysis_enabled and get_visible_url_names().get("upload")

    @cached_property
    def igv_links_enabled(self) -> bool:
        return any((self.analysis_enabled, settings.VARIANT_DETAILS_SHOW_ANNOTATION, settings.VARIANT_DETAILS_SHOW_SAMPLES))

    @cached_property
    def discordance_enabled(self) -> bool:
        return settings.DISCORDANCE_ENABLED


@timed_cache()
def get_settings_form_features() -> SettingsFormFeatures:
    return SettingsFormFeatures()


class SettingsOverrideForm(BaseModelForm):
    """ Warning: This hides some fields (in _hide_unused_fields) so not all fields from model will exist """
    class Meta:
        model = models.SettingsOverride
        fields = "__all__"
        widgets = {
            "email_weekly_updates": BlankNullBooleanSelect(),
            "email_discordance_updates": BlankNullBooleanSelect(),
            "variant_link_in_analysis_opens_new_tab": BlankNullBooleanSelect(),
            "tool_tips": BlankNullBooleanSelect(),
            "node_debug_tab": BlankNullBooleanSelect(),
            "import_messages": BlankNullBooleanSelect(),
            'default_sort_by_column': ModelSelect2(url='custom_column_autocomplete',
                                                  forward=['columns'],
                                                  attrs={'data-placeholder': 'Column...'}),
            'timezone': forms.Select(choices=[(None, "")] + [(tz, tz) for tz in settings.AVAILABLE_TZS], attrs={}),
            "allele_origin_exclude_filter": BlankNullBooleanSelect(),
            "grid_sample_label_template": TextInput(),
        }
        labels = {
            "email_weekly_updates": "Email Regular Updates",
            "email_discordance_updates": "Email Discordance Updates",
            "variant_link_in_analysis_opens_new_tab": "Variant Link in Analysis Opens New Tab",
            "tool_tips": "Tooltips",
            "tag_colors": "Tag Colours",
            "node_debug_tab": "Node Debug Tab",
            "import_messages": "Import Messages",
            "default_sort_by_column": "Default Sort by Column",
            "igv_port": "IGV Port",
            "default_genome_build": "Default Genome Build",
            "default_allele_origin": "Default Allele Origin Filter",
            "default_lab": "Default Lab",
            "timezone": "TimeZone (for downloads)",
            "allele_origin_focus": "Allele Origin focus",
            "allele_origin_exclude_filter": "Allele Origin (filter by default)",
            "grid_sample_label_template": "Grid Sample Label Template",
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['columns'].queryset = CustomColumnsCollection.filter_public()
        self.fields['default_genome_build'].queryset = GenomeBuild.builds_with_annotation()
        self._hide_unused_fields()

    def clean_grid_sample_label_template(self):
        data = self.cleaned_data["grid_sample_label_template"]
        if data:
            try:
                Sample._validate_sample_formatter_func(data)
            except (ValueError, KeyError) as e:
                raise ValidationError(e)
        return data

    def _hide_unused_fields(self):
        settings_config = get_settings_form_features()
        field_visibility = {
            # "email_weekly_updates": settings.DISCORDANCE_ENABLED, # this is also used for release note updates
            "email_discordance_updates": settings_config.discordance_enabled,
            "columns": settings_config.analysis_enabled,
            "default_sort_by_column": settings_config.analysis_enabled,
            "variant_link_in_analysis_opens_new_tab": settings_config.analysis_enabled,
            "tool_tips": settings_config.analysis_enabled,
            "node_debug_tab": settings_config.analysis_enabled,
            "tag_colors": settings_config.analysis_enabled,
            "import_messages": settings_config.upload_enabled,
            "igv_port": settings_config.igv_links_enabled,
            "grid_sample_label_template": settings_config.analysis_enabled,
        }

        for f, visible in field_visibility.items():
            if f in self.fields and not visible:
                del self.fields[f]


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
        fields = ['user', 'name', "genome_build"]
        widgets = {'user': HiddenInput(),
                   'name': TextInput()}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['genome_build'].choices = GenomeBuild.get_choices()


class CohortForm(forms.ModelForm):

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
                if ManualVariantEntry.get_entry_type(line) == ManualVariantEntryType.UNKNOWN:
                    msg = f"'{line}' does not match known patterns"
                    if " " in line:
                        msg += " (remember to only enter one variant per line)"
                    raise forms.ValidationError(msg)

        return cleaned_data


class UserCohortForm(forms.Form):
    cohort = forms.ModelChoiceField(queryset=Cohort.objects.all(),
                                    widget=ModelSelect2(url='cohort_autocomplete',
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

        genome_builds_qs = GenomeBuild.builds_with_annotation()
        build_field = self.fields['genome_build']
        build_field.queryset = genome_builds_qs

        self.genome_build_hidden = False

        if genome_builds_qs.count() < 2:
            self.genome_build_hidden = True
            build_field.widget = HiddenInput()


class GenomicIntervalsCollectionForm(forms.ModelForm, ROFormMixin):
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


class UserLabChoiceForm(forms.Form):
    lab = forms.ModelChoiceField(
        required=True,
        queryset=Lab.objects.all(),
        widget=ModelSelect2(url='lab_autocomplete',
                            attrs={'data-placeholder': 'Lab...'}))

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user")
        default_lab = kwargs.pop("default_lab")
        super().__init__(*args, **kwargs)
        self.fields['lab'].queryset = Lab.valid_labs_qs(user)
        if default_lab:
            lab_field = self.fields['lab']
            lab_field.initial = default_lab
            # lab_field.choices = [(default_lab.pk, default_lab)]
