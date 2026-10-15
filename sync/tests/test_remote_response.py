from unittest.mock import MagicMock, patch

import requests
from django.test import SimpleTestCase
from requests import Response

from sync.shariant.remote_response import (
    RETRY_ATTEMPTS,
    RemoteResponseError,
    request_json_with_retries,
)

PROXY_PAGE = b'<!DOCTYPE html>\n<!--[if lt IE 7]> <html class="no-js ie6 oldie" lang="en-US"> <![endif]-->'


def make_response(status_code: int, content_type: str, content: bytes, **headers) -> Response:
    response = Response()
    response.status_code = status_code
    response.headers.update({"Content-Type": content_type, **headers})
    response._content = content
    response.url = "https://shariant.example/classification/api/classifications/export"
    response.request = requests.Request("GET", response.url).prepare()
    return response


@patch("sync.shariant.remote_response.time.sleep")
class RequestJsonWithRetriesTest(SimpleTestCase):

    def test_proxy_html_page_retried_then_succeeds(self, _sleep):
        json_response = make_response(200, "application/json", b'{"records": []}')
        send = MagicMock(side_effect=[
            make_response(200, "text/html", PROXY_PAGE, Via="1.1 proxy"),
            requests.exceptions.ChunkedEncodingError(),
            json_response,
        ])
        self.assertIs(request_json_with_retries(send, "test"), json_response)

    def test_gives_up_naming_the_intermediary(self, _sleep):
        send = MagicMock(return_value=make_response(502, "text/html", PROXY_PAGE, Server="cloudflare", Via="1.1 proxy"))
        with self.assertRaises(RemoteResponseError) as cm:
            request_json_with_retries(send, "test")
        self.assertEqual(send.call_count, RETRY_ATTEMPTS)
        message = str(cm.exception)
        for expected in ("502", "Server=cloudflare", "Via=1.1 proxy", "<!DOCTYPE html>"):
            self.assertIn(expected, message)

    def test_client_error_not_retried(self, _sleep):
        send = MagicMock(return_value=make_response(401, "application/json", b'{"detail": "no"}'))
        with self.assertRaises(RemoteResponseError):
            request_json_with_retries(send, "test")
        self.assertEqual(send.call_count, 1)
