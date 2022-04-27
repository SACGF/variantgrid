from django import forms

from upload.models import UploadSettings, UploadedFileTypes


class UploadSettingsForm(forms.ModelForm):
    file_types = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
        choices=UploadedFileTypes.choices,
    )

    class Meta:
        model = UploadSettings
        exclude = ('user',)
        labels = {
            "time_filter_method": "",
            "time_filter_value": "",
        }

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        instance = kwargs.get("instance")
        initial = kwargs.get("initial") or {}
        kwargs["initial"] = initial

        usft_set = instance.uploadsettingsfiletype_set
        initial["file_types"] = list(usft_set.all().values_list("file_type", flat=True))
        super().__init__(*args, **kwargs)

    def save(self, commit=True):
        instance = super().save(commit=False)
        instance.user = self.user

        if commit:
            instance.save()
            usft_set = instance.uploadsettingsfiletype_set
            old_usft = set(usft_set.all().values_list("file_type", flat=True))
            new_usft = set(self.cleaned_data['file_types'])
            added = new_usft - old_usft
            removed = old_usft - new_usft
            for ft in added:
                usft_set.create(file_type=ft)
            if removed:
                usft_set.filter(file_type__in=removed).delete()

        return instance

    def clean_time_filter_value(self):
        data = self.cleaned_data['time_filter_value']
        if data < 1:
            raise forms.ValidationError("Minimum value of 1")
        return data
