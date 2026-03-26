from django.contrib.auth.models import User
from django.test import TestCase

from upload.forms import UploadSettingsForm
from upload.models import UploadSettings
from upload.models.models_enums import TimeFilterMethod


class TestUploadSettingsFormValidation(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="form_test_user", password="x")
        cls.upload_settings = UploadSettings.objects.create(user=cls.user)

    def _make_form(self, time_filter_value):
        return UploadSettingsForm(
            data={"time_filter_method": TimeFilterMethod.DAYS,
                  "time_filter_value": time_filter_value, "file_types": []},
            instance=self.upload_settings,
            user=self.user,
        )

    def test_time_filter_value_zero_is_invalid(self):
        form = self._make_form(0)
        self.assertFalse(form.is_valid())
        self.assertIn("time_filter_value", form.errors)

    def test_time_filter_value_one_is_valid(self):
        """Boundary: 1 is the minimum allowed value."""
        form = self._make_form(1)
        self.assertTrue(form.is_valid(), form.errors)
