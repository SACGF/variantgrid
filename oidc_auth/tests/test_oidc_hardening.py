from unittest.mock import Mock, patch

import jwt
from django.contrib.auth.models import User
from django.contrib.messages import get_messages
from django.contrib.messages.storage.fallback import FallbackStorage
from django.contrib.sessions.middleware import SessionMiddleware
from django.core.exceptions import PermissionDenied, SuspiciousOperation
from django.test import RequestFactory, TestCase, override_settings
from mozilla_django_oidc.contrib.drf import OIDCAuthentication
from mozilla_django_oidc.middleware import SessionRefresh
from requests import HTTPError, Response
from rest_framework.exceptions import AuthenticationFailed
from rest_framework.request import Request

from oidc_auth.backend import VariantGridOIDCAuthenticationBackend, provider_logout
from oidc_auth.oidc_error_handler import HandleOIDC400Middleware
from oidc_auth.session_refresh import VariantGridSessionRefresh

_OIDC_SETTINGS = {"OIDC_OP_TOKEN_ENDPOINT": "https://idp/token", "OIDC_OP_USER_ENDPOINT": "https://idp/userinfo",
                  "OIDC_RP_CLIENT_ID": "shariant", "OIDC_RP_CLIENT_SECRET": "secret",
                  "OIDC_REQUIRED_GROUP": "/variantgrid/shariant_production", "MAINTENANCE_MODE": False}

_PROD_CLAIMS = {"preferred_username": "oidc_user", "email": "oidc_user@example.com", "sub": "abc",
                "groups": ["/variantgrid/shariant_production", "/associations/org/lab"]}


def _http_error(status_code: int) -> HTTPError:
    response = Response()
    response.status_code = status_code
    return HTTPError("Get Token Error", response=response)


def _userinfo_response(claims: dict) -> Mock:
    response = Mock(headers={"content-type": "application/json"})
    response.json.return_value = claims
    return response


def _browser_login_backend():
    request = RequestFactory().get("/oidc/callback/")
    SessionMiddleware(lambda r: None).process_request(request)
    request._messages = FallbackStorage(request)
    backend = VariantGridOIDCAuthenticationBackend()
    backend.request = request
    return backend, request


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


@override_settings(**_OIDC_SETTINGS)
class OIDCBrowserLoginRefusalTest(TestCase):
    """ Browser callback refusals return an inactive user (mozilla's login_failure) and show a message """

    def test_email_is_escaped_in_html_message(self):
        backend, request = _browser_login_backend()
        email = "o'brien&co@example.com"
        user = User.objects.create(username="oidc_wrong_env")
        claims = {"preferred_username": "oidc_wrong_env", "email": email, "sub": "abc",
                  "groups": ["/variantgrid/shariant_test"]}
        backend.create_or_update(user, claims)

        message = str(next(iter(get_messages(request))))
        self.assertIn("<i>o&#x27;brien&amp;co@example.com</i>", message)
        self.assertIn('<br/>Please try out test environment <a href="https://test.shariant.org.au">', message)

    def test_refused_first_login_leaves_no_user(self):
        wrong_environment = {**_PROD_CLAIMS, "groups": ["/variantgrid/shariant_test"]}
        for claims, maintenance_mode in [(wrong_environment, False), (_PROD_CLAIMS, True)]:
            with self.subTest(maintenance_mode=maintenance_mode), override_settings(MAINTENANCE_MODE=maintenance_mode):
                backend, _ = _browser_login_backend()
                with patch("mozilla_django_oidc.auth.requests.get", return_value=_userinfo_response(claims)):
                    user = backend.get_or_create_user("access", "id", {})
                self.assertFalse(user.is_active)
                self.assertFalse(User.objects.filter(email=_PROD_CLAIMS["email"]).exists())

    def test_existing_user_without_labs_is_deactivated(self):
        backend, _ = _browser_login_backend()
        user = User.objects.create(username="oidc_user", email="oidc_user@example.com")
        with override_settings(OIDC_REQUIRED_GROUP=None):
            backend.create_or_update(user, {**_PROD_CLAIMS, "groups": ["/unrelated"]})
        user.refresh_from_db()
        self.assertFalse(user.is_active)


@override_settings(**_OIDC_SETTINGS, OIDC_API_EXTRA_CLIENT_IDS=["shariant-cli"])
class OIDCBearerTokenTest(TestCase):
    """ DRF's OIDCAuthentication path: no request on the backend, so refusals must be 401/403 rather than 500s """

    @staticmethod
    def _authenticate(claims, azp="shariant", aud="account"):
        access_token = jwt.encode({"azp": azp, "aud": aud}, "keycloak-signing-key-for-tests-only", algorithm="HS256")
        request = Request(RequestFactory().get("/classification/api/classifications/v3/record/",
                                               HTTP_AUTHORIZATION=f"Bearer {access_token}"))
        authentication = OIDCAuthentication(backend=VariantGridOIDCAuthenticationBackend())
        with patch("mozilla_django_oidc.auth.requests.get", return_value=_userinfo_response(claims)):
            return authentication.authenticate(request)

    def test_token_must_be_for_this_client(self):
        for azp, aud in [("shariant", "account"), ("shariant-cli", "account"), ("other", ["account", "shariant"])]:
            with self.subTest(azp=azp):
                user, _ = self._authenticate(_PROD_CLAIMS, azp=azp, aud=aud)
                self.assertTrue(user.is_active)

        with self.assertRaises(AuthenticationFailed):
            self._authenticate(_PROD_CLAIMS, azp="shariant-test")

    def test_refusals_are_403_without_creating_user(self):
        wrong_environment = {**_PROD_CLAIMS, "groups": ["/variantgrid/shariant_test", "/associations/org/lab"]}
        no_groups_claim = {k: v for k, v in _PROD_CLAIMS.items() if k != "groups"}
        for claims in [wrong_environment, no_groups_claim]:
            with self.subTest(groups=claims.get("groups")), self.assertRaises(PermissionDenied):
                self._authenticate(claims)

        with override_settings(MAINTENANCE_MODE=True), self.assertRaises(PermissionDenied):
            self._authenticate(_PROD_CLAIMS)
        self.assertFalse(User.objects.filter(email=_PROD_CLAIMS["email"]).exists())


class OIDCProviderLogoutTest(TestCase):

    @staticmethod
    def _request(session: dict):
        request = RequestFactory().post("/oidc/logout/")
        SessionMiddleware(lambda r: None).process_request(request)
        request.session.update(session)
        return request

    @override_settings(KEY_CLOAK_PROTOCOL_BASE="https://idp/protocol", LOGOUT_REDIRECT_URL="https://shariant")
    def test_password_login_skips_keycloak(self):
        self.assertEqual(provider_logout(self._request({})), "https://shariant")
        logout_url = provider_logout(self._request({"oidc_id_token": "token"}))
        self.assertTrue(logout_url.startswith("https://idp/protocol/logout?id_token_hint=token"))


@override_settings(OIDC_OP_AUTHORIZATION_ENDPOINT="https://idp/auth", OIDC_RP_CLIENT_ID="shariant")
class VariantGridSessionRefreshTest(TestCase):

    def test_only_api_prefixes_skip_refresh(self):
        middleware = VariantGridSessionRefresh(get_response=lambda request: None)
        factory = RequestFactory()
        for path in ["/api/v1/capabilities", "/classification/api/classifications/v3/record/"]:
            with self.subTest(path=path):
                self.assertFalse(middleware.is_refreshable_url(factory.get(path)))

        with patch.object(SessionRefresh, "is_refreshable_url", return_value=True):
            self.assertTrue(middleware.is_refreshable_url(factory.get("/variantopedia/view/api/")))
