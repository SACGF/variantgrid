""" create_manual_variant_entry_from_text: search links to it, so a GET only asks to confirm and a POST creates """
from unittest.mock import Mock, patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse


@override_settings(UPLOAD_ENABLED=True, VARIANT_MANUAL_CREATE=True, VARIANT_MANUAL_CREATE_BY_NON_ADMIN=True)
class CreateManualVariantFromTextTest(TestCase):

    def setUp(self):
        self.client.force_login(User.objects.create_user("create_manual_variant_user"))

    def _url(self, variants_text: str) -> str:
        return reverse("create_manual_variant_entry_from_text",
                       kwargs={"genome_build_name": "GRCh37", "variants_text": variants_text})

    @patch("annotation.views.create_manual_variants", return_value=Mock(pk=1))
    def test_get_confirms_post_creates(self, create_manual_variants):
        url = self._url("1:169519049 T>C")
        self.assertEqual(200, self.client.get(url).status_code)
        create_manual_variants.assert_not_called()

        self.assertEqual(302, self.client.post(url).status_code)
        create_manual_variants.assert_called_once()

    @patch("annotation.views.create_manual_variants")
    def test_one_variant_only(self, create_manual_variants):
        response = self.client.post(self._url("1:169519049 T>C\n1:169519050 A>G"))
        self.assertEqual(400, response.status_code)
        create_manual_variants.assert_not_called()
