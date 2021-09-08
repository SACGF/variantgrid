from django import forms

from upload.models import UploadSettings


class UploadSettingsForm(forms.ModelForm):

    class Meta:
        model = UploadSettings
        exclude = ('user',)

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user")
        super().__init__(*args, **kwargs)

    def save(self, commit=True):
        instance = super().save(commit=False)
        instance.user = self.user

        if commit:
            instance.save()

        return instance

    def clean_time_filter_value(self):
        data = self.cleaned_data['time_filter_value']
        if data < 1:
            raise forms.ValidationError("Minimum value of 1")
        return data
