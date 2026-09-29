from django.conf import settings
from django.test import TestCase
from django.urls import reverse

from library.django_utils.login_required_middleware import login_not_required_include


class SiteLoginRequiredMiddlewareTest(TestCase):

    def test_drf_view_is_left_to_drf(self):
        """ DRF marks its views login_not_required: DRF answers, rather than a redirect to the login page """
        self.assertIn(self.client.get(reverse("classification_api")).status_code, (401, 403))

    def test_plain_view_under_an_api_prefix_needs_login(self):
        """ Only the view decides - a non-DRF view under /classification/api/ is not public for its path """
        url = reverse("imported_allele_info_datatables")
        response = self.client.get(url)
        self.assertRedirects(response, f"{settings.LOGIN_URL}?next={url}", fetch_redirect_response=False)

    def test_admin_goes_to_the_site_login(self):
        """ Django would send these to the admin's password form: its views name it as their login_url, and the
            admin login itself is login_not_required """
        for url in (reverse("admin:index"), reverse("admin:login")):
            with self.subTest(url=url):
                response = self.client.get(url)
                self.assertRedirects(response, f"{settings.LOGIN_URL}?next={url}", fetch_redirect_response=False)

    def test_login_not_required_include_marks_nested_urlconfs(self):
        """ registration.backends.default.urls includes registration.auth_urls; OIDC is only mounted on Shariant """
        for urlconf_module in ("registration.backends.default.urls", "mozilla_django_oidc.urls"):
            urlconf, _app_name, _namespace = login_not_required_include(urlconf_module)
            with self.subTest(urlconf=urlconf_module):
                self.assertTrue(self._all_callbacks(urlconf.urlpatterns))
                for callback in self._all_callbacks(urlconf.urlpatterns):
                    self.assertFalse(callback.login_required)

    def _all_callbacks(self, url_patterns) -> list:
        callbacks = []
        for entry in url_patterns:
            if hasattr(entry, "url_patterns"):
                callbacks.extend(self._all_callbacks(entry.url_patterns))
            else:
                callbacks.append(entry.callback)
        return callbacks
