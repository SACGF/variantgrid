from django.contrib.auth.models import User
from django.contrib.messages import get_messages
from django.contrib.messages.storage.fallback import FallbackStorage
from django.contrib.sessions.middleware import SessionMiddleware
from django.core.exceptions import SuspiciousOperation
from django.test import RequestFactory, TestCase, override_settings
from django.urls import ResolverMatch
from requests import HTTPError, Response

from oidc_auth.backend import VariantGridOIDCAuthenticationBackend, provider_logout
from oidc_auth.oidc_error_handler import HandleOIDC400Middleware
from snpdb.models.models_user_settings import UserSettingsOverride

_OIDC_CLIENT_SETTINGS = {"OIDC_OP_TOKEN_ENDPOINT": "https://idp/token", "OIDC_OP_USER_ENDPOINT": "https://idp/userinfo",
                         "OIDC_RP_CLIENT_ID": "variantgrid", "OIDC_RP_CLIENT_SECRET": "secret"}


def _http_error(status_code: int) -> HTTPError:
    response = Response()
    response.status_code = status_code
    return HTTPError("Get Token Error", response=response)


class HandleOIDC400MiddlewareTest(TestCase):
    """ Only a stale/replayed /oidc/ callback is swallowed into a redirect """

    def setUp(self):
        self.middleware = HandleOIDC400Middleware(get_response=lambda request: None)
        self.factory = RequestFactory()

    def _callback_view(self, oidc_states=None, url_name="oidc_authentication_callback"):
        """ process_view on a callback for state "abc", with the session's oidc_states (None = no key) """
        request = self.factory.get("/oidc/callback/", {"state": "abc", "code": "abc"})
        request.resolver_match = ResolverMatch(lambda r: None, (), {}, url_name=url_name)
        SessionMiddleware(lambda r: None).process_request(request)
        if oidc_states is not None:
            request.session["oidc_states"] = oidc_states
        return self.middleware.process_view(request, None, (), {})

    def test_unknown_callback_state_redirects_to_root(self):
        response = self._callback_view(oidc_states={"other": {}})
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, "/")

    def test_callback_left_to_view(self):
        """ A known state, no oidc_states at all (the view's login_failure) or another URL isn't intercepted """
        self.assertIsNone(self._callback_view(oidc_states={"abc": {}}))
        self.assertIsNone(self._callback_view())
        self.assertIsNone(self._callback_view(oidc_states={"other": {}}, url_name="index"))

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


@override_settings(**_OIDC_CLIENT_SETTINGS, OIDC_REQUIRED_GROUP="/variantgrid/shariant_production")
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


@override_settings(**_OIDC_CLIENT_SETTINGS)
class OIDCUserLookupTest(TestCase):
    """ filter_users_by_claims: the stored Keycloak sub wins, unlinked users are found by username or verified email,
        and a user linked to another sub is never taken over """

    def setUp(self):
        self.backend = VariantGridOIDCAuthenticationBackend()
        self.alice = User.objects.create(username="alice", email="alice@example.com")
        self.bob = User.objects.create(username="bob", email="bob@example.com")

    @staticmethod
    def _link(user, sub):
        UserSettingsOverride.objects.update_or_create(user=user, defaults={"oauth_sub": sub})

    def _lookup(self, **claims):
        return list(self.backend.filter_users_by_claims(claims))

    def test_sub_match_wins(self):
        self._link(self.alice, "sub-alice")
        # Renamed in Keycloak: still alice's account, not bob's
        self.assertEqual(self._lookup(sub="sub-alice", preferred_username="bob", email="bob@example.com",
                                      email_verified=True), [self.alice])

    def test_unlinked_user_found_by_username(self):
        """ Existing users with no stored sub still log in (create_or_update then stores it) """
        self.assertEqual(self._lookup(sub="sub-new", preferred_username="alice"), [self.alice])

    def test_email_match_needs_verified_email(self):
        claims = {"sub": "sub-new", "preferred_username": "someone_else", "email": "ALICE@example.com"}
        self.assertEqual(self._lookup(**claims), [])
        self.assertEqual(self._lookup(**claims, email_verified=False), [])
        self.assertEqual(self._lookup(**claims, email_verified=True), [self.alice])

    def test_user_linked_to_another_sub_is_refused(self):
        self._link(self.alice, "sub-alice")
        with self.assertRaises(SuspiciousOperation):
            self._lookup(sub="sub-other", preferred_username="alice")
        with self.assertRaises(SuspiciousOperation):
            self._lookup(sub="sub-other", preferred_username="nobody", email="alice@example.com", email_verified=True)


@override_settings(**_OIDC_CLIENT_SETTINGS, OIDC_REQUIRED_GROUP="/variantgrid/shariant_production")
class OIDCMissingClaimsTest(TestCase):

    def test_missing_groups_claim_deactivates(self):
        request = RequestFactory().get("/oidc/callback/")
        SessionMiddleware(lambda r: None).process_request(request)
        request._messages = FallbackStorage(request)
        backend = VariantGridOIDCAuthenticationBackend()
        backend.request = request

        user = User.objects.create(username="oidc_no_groups")
        claims = {"preferred_username": "oidc_no_groups", "email": "nogroups@example.com", "sub": "abc"}
        backend.create_or_update(user, claims)
        self.assertFalse(user.is_active)

    @override_settings(KEY_CLOAK_PROTOCOL_BASE="https://idp/protocol", LOGOUT_REDIRECT_URL="https://vg/")
    def test_logout_without_id_token(self):
        """ A password (ModelBackend) login has no oidc_id_token in its session """
        request = RequestFactory().post("/oidc_logout")
        SessionMiddleware(lambda r: None).process_request(request)
        self.assertEqual(provider_logout(request), "https://idp/protocol/logout")
