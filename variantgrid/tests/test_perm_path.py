from collections import defaultdict
from unittest import mock

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.http import HttpResponse
from django.test import RequestFactory, SimpleTestCase, TestCase, override_settings
from rest_framework import routers

from classification.views.classification_view import ClassificationView
from patients.views_rest import PatientViewSet
from variantgrid.perm_path import deprecated_path, get_visible_url_names, router_urls

REGISTER = defaultdict(lambda: True, {"api_patient-list": False})


class UrlNameSettingsTest(TestCase):
    """ lab_members_tab follows LAB_HEAD_MANAGE_MEMBERS when the register is read, not when settings loaded """

    def setUp(self):
        get_visible_url_names.cache_clear()
        self.addCleanup(get_visible_url_names.cache_clear)

    @override_settings(LAB_HEAD_MANAGE_MEMBERS=False)
    def test_setting_off_hides_url(self):
        self.assertFalse(get_visible_url_names()["lab_members_tab"])

    @override_settings(LAB_HEAD_MANAGE_MEMBERS=True)
    def test_setting_on_shows_url(self):
        self.assertTrue(get_visible_url_names()["lab_members_tab"])


class RouterUrlsTest(TestCase):
    """ The register is read when a URLconf is imported, so build a router here rather than reversing
        through the global one """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user("perm_path_user")
        cls.super_user = User.objects.create_superuser("perm_path_super_user")

    @staticmethod
    def _url_patterns() -> dict:
        router = routers.DefaultRouter()
        router.register(r'api/v1/patient', PatientViewSet, basename='api_patient')
        return {url.name: url for url in router_urls(router)}

    def _get(self, url, user, **kwargs):
        request = RequestFactory().get("/patients/api/v1/patient/")
        request.user = user
        return url.callback(request, **kwargs)

    @override_settings(URLS_NAME_REGISTER=REGISTER)
    def test_unregistered_name_requires_superuser(self):
        url = self._url_patterns()["api_patient-list"]
        with self.assertRaises(PermissionDenied):
            self._get(url, self.user)
        self.assertEqual(self._get(url, self.super_user).status_code, 200)

    @override_settings(URLS_NAME_REGISTER=REGISTER)
    def test_registered_name_passes_through(self):
        """ api_patient-detail is left True, so the viewset handles a non-superuser itself """
        url = self._url_patterns()["api_patient-detail"]
        self.assertEqual(self._get(url, self.user, pk=0).status_code, 404)


class DeprecatedPathTest(SimpleTestCase):

    def test_reports_and_passes_through(self):
        url = deprecated_path("old/<int:pk>", lambda request, pk: HttpResponse(str(pk)), name="old_view")
        request = RequestFactory().get("/old/5?x=1")
        with mock.patch("variantgrid.perm_path.report_message") as report_message:
            response = url.callback(request, pk=5)
        self.assertEqual(response.content, b"5")
        report_message.assert_called_once()
        self.assertIn("old_view", report_message.call_args.args[0])
        self.assertEqual(report_message.call_args.kwargs["extra_data"]["target"], "/old/5?x=1")

    def test_keeps_api_view_exemptions(self):
        """ A deprecated DRF view stays callable by external clients without a CSRF token or session """
        url = deprecated_path("api/v2/", ClassificationView.as_view(api_version=2), name="old_api")
        self.assertTrue(url.callback.csrf_exempt)
        self.assertFalse(url.callback.login_required)
