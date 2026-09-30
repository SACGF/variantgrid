from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from library.django_utils.view_utils import login_not_required_include


class LoginRequiredTest(TestCase):

    def test_drf_view_is_left_to_drf(self):
        """ DRF marks its views login_not_required: DRF answers, rather than a redirect to the login page """
        self.assertIn(self.client.get(reverse("classification_api")).status_code, (401, 403))

    def test_plain_view_under_an_api_prefix_needs_login(self):
        """ Only the view decides - a non-DRF view under /classification/api/ is not public for its path """
        url = reverse("imported_allele_info_datatables")
        response = self.client.get(url)
        self.assertRedirects(response, f"{settings.LOGIN_URL}?next={url}", fetch_redirect_response=False)

    def test_admin_login_goes_to_the_site_login(self):
        """ Django marks the admin login login_not_required, and admin views name it as their login_url """
        admin_login = reverse("admin:login")
        response = self.client.get(admin_login)
        self.assertRedirects(response, f"{settings.LOGIN_URL}?next={admin_login}", fetch_redirect_response=False)

        response = self.client.get(reverse("admin:index"))
        self.assertTrue(response.url.startswith(admin_login))

    def test_staff_only(self):
        """ Anonymous users go to the site login, a logged-in non-staff user to the staff_only page """
        url = reverse("view_sequencer", kwargs={"pk": 1})
        response = self.client.get(url)
        self.assertRedirects(response, f"{settings.LOGIN_URL}?next={url}", fetch_redirect_response=False)

        self.client.force_login(User.objects.create_user("not_staff"))
        response = self.client.get(url)
        self.assertRedirects(response, reverse("staff_only"), fetch_redirect_response=False)

    def test_login_not_required_include_marks_nested_urlconfs(self):
        """ registration.backends.default.urls includes registration.auth_urls; OIDC is only mounted on Shariant """
        for urlconf_module in ("registration.backends.default.urls", "mozilla_django_oidc.urls"):
            urlconf, _app_name, _namespace = login_not_required_include(urlconf_module)
            with self.subTest(urlconf=urlconf_module):
                callbacks = self._all_callbacks(urlconf.urlpatterns)
                self.assertTrue(callbacks)
                for callback in callbacks:
                    self.assertFalse(callback.login_required)

    def _all_callbacks(self, url_patterns) -> list:
        callbacks = []
        for entry in url_patterns:
            if hasattr(entry, "url_patterns"):
                callbacks.extend(self._all_callbacks(entry.url_patterns))
            else:
                callbacks.append(entry.callback)
        return callbacks
