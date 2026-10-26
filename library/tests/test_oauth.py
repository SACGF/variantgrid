from django.test import SimpleTestCase

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
