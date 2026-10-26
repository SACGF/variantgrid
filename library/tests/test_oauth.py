import time
from unittest.mock import patch

from django.test import SimpleTestCase
from oauthlib.oauth2 import InvalidGrantError
from requests import Request
from requests_oauthlib.oauth2_session import OAuth2Session

from library.oauth import ServerAuth


class ServerAuthUrlTest(SimpleTestCase):

    @staticmethod
    def _server_auth(host: str) -> ServerAuth:
        return ServerAuth(username="user", password="password", host=host)

    def test_url_keeps_host_sub_path(self):
        server_auth = self._server_auth("https://example.com/variantgrid")
        expected = "https://example.com/variantgrid/classification/api/x"
        self.assertEqual(server_auth.url("/classification/api/x"), expected)
        self.assertEqual(server_auth.url("classification/api/x"), expected)

    def test_non_http_host_rejected(self):
        with self.assertRaises(ValueError):
            self._server_auth("file:///etc/passwd").url("x")

    def test_http_host_needs_allow_http(self):
        with self.assertRaises(ValueError):
            self._server_auth("http://example.com").url("x")
        allowed = ServerAuth(username="user", password="password", host="http://example.com", allow_http=True)
        self.assertEqual(allowed.url("x"), "http://example.com/x")
        self.assertEqual(self._server_auth("http://localhost:8080").url("x"), "http://localhost:8080/x")


@patch("library.oauth.requests.get")
class ServerAuthTokenTest(SimpleTestCase):
    """ Token requests are mocked at OAuth2Session, so no Keycloak is contacted """

    def setUp(self):
        self.tokens_issued = 0
        fetch_patcher = patch.object(OAuth2Session, "fetch_token", autospec=True, side_effect=self._issue_token)
        refresh_patcher = patch.object(OAuth2Session, "refresh_token", autospec=True, side_effect=self._issue_token)
        self.fetch_token = fetch_patcher.start()
        self.refresh_token = refresh_patcher.start()
        self.addCleanup(fetch_patcher.stop)
        self.addCleanup(refresh_patcher.stop)
        self.server_auth = ServerAuth(username="user", password="password", host="https://shariant.example.com",
                                      oauth_url="https://keycloak.example.com/token", client_id="client")

    def _issue_token(self, session: OAuth2Session, *args, **kwargs):
        self.tokens_issued += 1
        session.token = {"access_token": f"token{self.tokens_issued}", "token_type": "Bearer",
                         "refresh_token": "refresh", "expires_at": time.time() + 300}
        return session.token

    def _sent_access_token(self, mock_get) -> str:
        request = mock_get.call_args.kwargs["auth"](Request("GET", "https://shariant.example.com").prepare())
        return request.headers["Authorization"].removeprefix("Bearer ")

    def _get_after_expiry(self):
        self.server_auth.get("classification/api/x")
        with patch("library.oauth.time.time", return_value=time.time() + 600):
            self.server_auth.get("classification/api/x")

    def test_one_login_for_several_requests(self, mock_get):
        for _ in range(3):
            self.server_auth.get("classification/api/x")
        self.fetch_token.assert_called_once()
        self.refresh_token.assert_not_called()
        self.assertEqual(self._sent_access_token(mock_get), "token1")

    def test_expired_token_refreshed(self, mock_get):
        self._get_after_expiry()
        self.fetch_token.assert_called_once()
        self.refresh_token.assert_called_once()
        self.assertEqual(self._sent_access_token(mock_get), "token2")

    def test_refused_refresh_logs_in_again(self, mock_get):
        self.refresh_token.side_effect = InvalidGrantError()
        self._get_after_expiry()
        self.assertEqual(self.fetch_token.call_count, 2)
        self.assertEqual(self._sent_access_token(mock_get), "token2")
