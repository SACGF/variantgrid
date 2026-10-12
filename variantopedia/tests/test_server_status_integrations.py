from contextlib import contextmanager

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse

from annotation.fake_data import create_fake_variants, get_fake_annotation_version
from annotation.models import AnnotationRangeLock, AnnotationRun
from library.django_utils.unittest_utils import URLTestCase
from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    integration_status_signal,
)
from library.tests.test_integration_status import temporarily_connected
from snpdb.models import GenomeBuild, Variant
from variantopedia.views_server_status import _highest_variant_annotation_status


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


class ServerStatusAnnotationCheckTest(TestCase):
    """ The highest-variant annotation check on Server Status - see _highest_variant_annotation_status """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch37()
        cls.vav = get_fake_annotation_version(cls.genome_build).variant_annotation_version
        create_fake_variants(cls.genome_build)

    def test_unannotated_with_no_run_is_danger(self):
        status = _highest_variant_annotation_status()
        self.assertEqual("danger", status["status"])
        self.assertEqual("Not annotated, no AnnotationRun!", status["message"])

    def test_unannotated_with_a_covering_run_is_warning(self):
        variants = Variant.objects.order_by("pk")
        range_lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=variants.first(),
                                                        max_variant=variants.last())
        AnnotationRun.objects.create(annotation_range_lock=range_lock)
        status = _highest_variant_annotation_status()
        self.assertEqual("warning", status["status"])
        self.assertTrue(status["message"].startswith("AnnotationRuns: "))


@override_settings(CELERY_ENABLED=False)
class ServerStatusPostTest(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.admin_user = User.objects.create_superuser("server_status_post_admin")

    def test_action_redirects_so_a_refresh_does_not_repeat_it(self):
        self.client.force_login(self.admin_user)
        response = self.client.post(reverse('server_status'), {"action": "Test Message Branding"})
        self.assertRedirects(response, reverse('server_status'), fetch_redirect_response=False)
