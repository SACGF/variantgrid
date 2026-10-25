from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse


class WaitForTaskRedirectTest(TestCase):
    """ wait_for_task only redirects within the site """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="wait_for_task_user")

    def _get(self, redirect_url):
        self.client.force_login(self.user)
        kwargs = {"celery_task": "no-such-task", "sleep_ms": 1000, "redirect_url": redirect_url}
        return self.client.get(reverse("wait_for_task", kwargs=kwargs))

    def test_other_site_refused(self):
        self.assertEqual(self._get("https://example.com/").status_code, 400)

    def test_same_site_allowed(self):
        self.assertEqual(self._get("/snpdb/").status_code, 200)
