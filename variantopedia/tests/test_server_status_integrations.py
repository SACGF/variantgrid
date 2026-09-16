from contextlib import contextmanager

from django.contrib.auth.models import User
from django.test import override_settings
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase
from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    integration_status_signal,
)
from library.tests.test_integration_status import temporarily_connected


@contextmanager
def no_providers_registered():
    """ A deployment with nothing registered - the section should vanish entirely """
    registered = integration_status_signal.receivers
    integration_status_signal.receivers = []
    integration_status_signal.sender_receivers_cache.clear()
    try:
        yield
    finally:
        integration_status_signal.receivers = registered
        integration_status_signal.sender_receivers_cache.clear()


@override_settings(CELERY_ENABLED=False)
class ServerStatusIntegrationsTest(URLTestCase):
    """ The Integrations section renders inline on Server Status - see library/integration_status.py """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.admin_user = User.objects.create_superuser(f"test_user_{__file__}_admin")

    def _server_status_html(self) -> str:
        self.client.force_login(self.admin_user)
        response = self.client.get(reverse('server_status'))
        self.assertEqual(200, response.status_code)
        return response.content.decode()

    def test_renders_an_integration_that_has_never_run(self):
        def provider(sender, **kwargs):
            return IntegrationStatus(name="Never Run Integration",
                                     details=[IntegrationDetail(label="Last Run")])

        with temporarily_connected(provider):
            html = self._server_status_html()

        self.assertIn("Integrations", html)
        self.assertIn("Never Run Integration", html)
        self.assertIn("Never run", html)

    def test_section_hidden_when_nothing_registered(self):
        with no_providers_registered():
            html = self._server_status_html()
        self.assertNotIn("Integrations", html)
