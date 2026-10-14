from django.conf import settings
from django.contrib.sites.models import Site
from django.test import SimpleTestCase, TestCase, override_settings

from snpdb.checks import OIDC_BACKEND
from variantgrid.deployment_validation.settings_check import (
    check_settings_references,
    check_system_checks,
)


class SystemChecksAsDeploymentCheckTest(SimpleTestCase):

    def test_current_settings_pass(self):
        for name, result in check_system_checks().items():
            with self.subTest(name):
                self.assertTrue(result["valid"], result)

    @override_settings(AUTHENTICATION_BACKENDS=(*settings.AUTHENTICATION_BACKENDS, OIDC_BACKEND),
                       DISCORDANCE_EMAIL="updates@example.org", SEND_EMAILS=False)
    def test_error_is_invalid_and_warning_is_a_warning(self):
        results = check_system_checks()
        self.assertFalse(results["variantgrid.E005"]["valid"])
        self.assertIn("USE_OIDC is False", results["variantgrid.E005"]["fix"])
        self.assertTrue(results["variantgrid.W001"]["valid"])
        self.assertIn("SEND_EMAILS is False", results["variantgrid.W001"]["warning"])


class SettingsReferencesCheckTest(TestCase):

    def test_rows_named_by_settings_must_exist(self):
        site, _ = Site.objects.get_or_create(pk=settings.SITE_ID, defaults={"domain": "example.org", "name": "example"})
        with override_settings(LAB_EXTERNAL_MEMBERSHIP={"no_such_org/no_such_lab": "Ask IT"}):
            results = check_settings_references()
        self.assertTrue(results["SITE_ID"]["valid"])
        self.assertFalse(results["LAB_EXTERNAL_MEMBERSHIP['no_such_org/no_such_lab']"]["valid"])

        with override_settings(SITE_ID=site.pk + 1000):
            self.assertFalse(check_settings_references()["SITE_ID"]["valid"])
