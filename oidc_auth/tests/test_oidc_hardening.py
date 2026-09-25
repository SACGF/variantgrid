from django.contrib.auth.models import User
from django.contrib.messages import get_messages
from django.contrib.messages.storage.fallback import FallbackStorage
from django.contrib.sessions.middleware import SessionMiddleware
from django.core.exceptions import SuspiciousOperation
from django.test import RequestFactory, TestCase, override_settings
from requests import HTTPError, Response

from oidc_auth.backend import VariantGridOIDCAuthenticationBackend
from oidc_auth.oidc_error_handler import HandleOIDC400Middleware


def _http_error(status_code: int) -> HTTPError:
    response = Response()
    response.status_code = status_code
    return HTTPError("Get Token Error", response=response)


class HandleOIDC400MiddlewareTest(TestCase):
    """ Only the provider rejecting a stale code (400) on /oidc/ is swallowed into a redirect """

    def setUp(self):
        self.middleware = HandleOIDC400Middleware(get_response=lambda request: None)
        self.factory = RequestFactory()

    def test_rejected_token_on_oidc_path_redirects_to_root(self):
        request = self.factory.get("/oidc/callback/")
        response = self.middleware.process_exception(request, _http_error(400))
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, "/")

    def test_other_exceptions_are_left_to_django(self):
        request = self.factory.get("/oidc/callback/")
        for exception in [SuspiciousOperation("state mismatch"), _http_error(500), ValueError("bug")]:
            with self.subTest(exception=exception):
                self.assertIsNone(self.middleware.process_exception(request, exception))

        self.assertIsNone(self.middleware.process_exception(self.factory.get("/snpdb/"), _http_error(400)))


@override_settings(OIDC_OP_TOKEN_ENDPOINT="https://idp/token", OIDC_OP_USER_ENDPOINT="https://idp/userinfo",
                   OIDC_RP_CLIENT_ID="variantgrid", OIDC_RP_CLIENT_SECRET="secret",
                   OIDC_REQUIRED_GROUP="/variantgrid/shariant_production")
class OIDCWrongEnvironmentMessageTest(TestCase):

    def test_email_is_escaped_in_html_message(self):
        request = RequestFactory().get("/oidc/callback/")
        SessionMiddleware(lambda r: None).process_request(request)
        request._messages = FallbackStorage(request)
        backend = VariantGridOIDCAuthenticationBackend()
        backend.request = request

        email = "o'brien&co@example.com"
        user = User.objects.create(username="oidc_wrong_env")
        claims = {"preferred_username": "oidc_wrong_env", "email": email, "sub": "abc",
                  "groups": ["/variantgrid/shariant_test"]}
        backend.create_or_update(user, claims)

        message = str(next(iter(get_messages(request))))
        self.assertIn("<i>o&#x27;brien&amp;co@example.com</i>", message)
        self.assertIn('<br/>Please try out test environment <a href="https://test.shariant.org.au">', message)
