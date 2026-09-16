from collections import defaultdict

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import RequestFactory, TestCase, override_settings
from rest_framework import routers

from patients.views_rest import PatientViewSet
from variantgrid.perm_path import router_urls

REGISTER = defaultdict(lambda: True, {"api_patient-list": False})


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
