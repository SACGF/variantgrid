"""
sapath#443 - the capabilities endpoint a client asks before choosing which calls to make
"""
from django.contrib.auth.models import User
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase


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
        self.assertIn("version", data)
        self.assertIn("patients", data["features"])
        self.assertIn("dragen_tso500_combined_variant_output", data["upload_file_types"])
        self.assertNotIn("liftover", data["upload_file_types"])

    def test_anonymous_is_refused(self):
        response = self.client.get(reverse("api_capabilities"))
        self.assertEqual(response.status_code, status.HTTP_401_UNAUTHORIZED)
