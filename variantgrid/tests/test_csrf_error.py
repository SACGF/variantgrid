from django.contrib.auth.models import User
from django.test import Client, TestCase


class CsrfErrorTest(TestCase):
    def test_csrf_failure_is_forbidden_not_server_error(self):
        user = User.objects.create_user("csrf_error_user")
        client = Client(enforce_csrf_checks=True)
        client.force_login(user)
        response = client.post("/", {})
        self.assertEqual(response.status_code, 403)
        self.assertIn(b"Your session token has changed", response.content)
