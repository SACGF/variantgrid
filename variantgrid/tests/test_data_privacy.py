"""
#2852 - the public Data & privacy page, shown (and linked from the login page) only where settings.DATA_PRIVACY is set
"""
from django.test import TestCase, override_settings
from django.urls import reverse

DATA_PRIVACY = {"hosting": "agent_test hosting region", "outside_eu": True}


class DataPrivacyPageTest(TestCase):

    @override_settings(DATA_PRIVACY=DATA_PRIVACY, USER_CREATE_ORG_LABS={"agent_test_org": "user_group_%(username)s"})
    def test_configured_page_is_public_and_linked_from_login(self):
        response = self.client.get(reverse("data_privacy"))
        self.assertContains(response, "agent_test hosting region")
        self.assertContains(response, "GDPR")
        self.assertContains(response, "Getting a lab for your institution")
        self.assertContains(self.client.get(reverse("auth_login")), reverse("data_privacy"))

    @override_settings(DATA_PRIVACY=None)
    def test_unconfigured_page_404s_and_is_not_linked(self):
        self.assertEqual(self.client.get(reverse("data_privacy")).status_code, 404)
        self.assertNotContains(self.client.get(reverse("auth_login")), reverse("data_privacy"))
