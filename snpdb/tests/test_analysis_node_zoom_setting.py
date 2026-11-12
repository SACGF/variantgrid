from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from snpdb.forms import UserSettingsOverrideForm
from snpdb.models import Country, Lab, Organization
from snpdb.models.models_user_settings import LabUserSettingsOverride, UserSettingsOverride


class SetAnalysisNodeZoomTest(TestCase):
    """ The analysis settings zoom slider saves to the user's own override, and clearing it falls back
        to the lab's value (#2784) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        org = Organization.objects.create(name="ZoomOrg", group_name="zoomorg")
        country = Country.objects.get_or_create(name="CountryA")[0]
        cls.lab = Lab.objects.create(name="ZoomLab", organization=org, city="CityA",
                                     country=country, group_name="zoomorg/zoomlab")
        cls.user = User.objects.create(username="node_zoom_user")
        cls.user.groups.add(cls.lab.group)
        UserSettingsOverride.objects.update_or_create(user=cls.user, defaults={"default_lab": cls.lab})
        LabUserSettingsOverride.objects.update_or_create(lab=cls.lab, defaults={"analysis_node_zoom": 150})

    def setUp(self):
        self.client.force_login(self.user)

    def _post(self, zoom):
        return self.client.post(reverse("set_analysis_node_zoom"), {"zoom": zoom})

    def _user_zoom(self):
        return UserSettingsOverride.objects.get(user=self.user).analysis_node_zoom

    def test_set_then_clear_falls_back_to_lab(self):
        response = self._post("120")
        self.assertEqual(response.json(), {"zoom": 120})
        self.assertEqual(self._user_zoom(), 120)

        response = self._post("")
        self.assertEqual(response.json(), {"zoom": 150})
        self.assertIsNone(self._user_zoom())

    def test_rejects_zoom_off_the_slider_steps(self):
        for zoom in ["1000", "125", "abc"]:
            response = self._post(zoom)
            self.assertEqual(response.status_code, 400, zoom)
        self.assertIsNone(self._user_zoom())


class AnalysisNodeZoomSettingsFormTest(TestCase):
    """ The settings form's slider sits behind a checkbox: unticked saves null whatever the slider posts """

    def test_checkbox_sets_or_clears(self):
        user = User.objects.create(username="node_zoom_form_user")
        override = UserSettingsOverride.objects.update_or_create(user=user, defaults={"analysis_node_zoom": 150})[0]

        form = UserSettingsOverrideForm({"analysis_node_zoom": "120"}, instance=override)
        self.assertTrue(form.is_valid(), form.errors)
        self.assertIsNone(form.cleaned_data["analysis_node_zoom"])

        form = UserSettingsOverrideForm({"analysis_node_zoom-set": "on", "analysis_node_zoom": "120"}, instance=override)
        self.assertTrue(form.is_valid(), form.errors)
        self.assertEqual(form.cleaned_data["analysis_node_zoom"], 120)
