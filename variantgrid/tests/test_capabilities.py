"""
sapath#443 - the capabilities endpoint (/api/v1/capabilities) a client asks before choosing which calls to make.
Public, with the exact commit only for an authenticated caller (#2732)
"""
from django.contrib.auth.models import User
from django.test import SimpleTestCase, override_settings
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from variantgrid.views_rest import release_version


@override_settings(VARIANTGRID_VERSION="vg4.0-12-gc174556")
class CapabilitiesAPITest(APITestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="capabilities_user")

    def test_capabilities(self):
        self.client.force_authenticate(user=self.user)
        response = self.client.get(reverse("api_capabilities"))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertEqual(data["version"], "vg4.0-12-gc174556")
        self.assertTrue(data["git_hash"])
        self.assertIn("patients", data["features"])
        self.assertIn("dragen_tso500_combined_variant_output", data["upload_file_types"])
        self.assertIn("dragen_tso500_metrics_output", data["upload_file_types"])
        self.assertNotIn("liftover", data["upload_file_types"])

    def test_anonymous_gets_capabilities_without_commit(self):
        response = self.client.get(reverse("api_capabilities"))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = response.json()
        self.assertEqual(data["version"], "vg4.0")
        self.assertNotIn("git_hash", data)
        self.assertIn("patients", data["features"])
        self.assertIn("vcf", data["upload_file_types"])


class ReleaseVersionTest(SimpleTestCase):

    def test_release_version(self):
        self.assertEqual(release_version("vg4.0-12-gc174556"), "vg4.0")
        self.assertEqual(release_version("vg4-gc174556"), "vg4")
        self.assertEqual(release_version("vg4.0"), "vg4.0")
