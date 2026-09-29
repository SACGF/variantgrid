from django.conf import settings
from django.test import TestCase
from django.urls import reverse


class PublicPathsLoginRequiredMiddlewareTest(TestCase):

    def test_public_path_is_left_to_drf(self):
        """ /classification/api/ is in PUBLIC_PATHS: DRF answers, rather than a redirect to the login page """
        self.assertIn(self.client.get(reverse("classification_api")).status_code, (401, 403))

    def test_admin_goes_to_the_site_login(self):
        """ Django would send these to the admin's password form: its views name it as their login_url, and the
            admin login itself is login_not_required """
        for url in (reverse("admin:index"), reverse("admin:login")):
            with self.subTest(url=url):
                response = self.client.get(url)
                self.assertRedirects(response, f"{settings.LOGIN_URL}?next={url}", fetch_redirect_response=False)
