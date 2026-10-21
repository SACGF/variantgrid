"""
Requests to a remote VariantGrid (Shariant) for the sync runners: `request_json_with_retries` sends one, retries what
an intermediary can cause - a dropped connection, a timeout, a truncated body, a 5xx or an HTML page in place of
JSON - and raises `RemoteResponseError` naming the status, proxy headers and start of the body when it gives up (#1872).
"""
import logging
import time
from collections.abc import Callable

import requests
from requests import Response

RETRY_ATTEMPTS = 3
RETRY_DELAY_SECS = 30

# headers that name whatever answered in place of the remote - a proxy, firewall or CDN
DIAGNOSTIC_HEADERS = ('Server', 'Via', 'cf-ray', 'X-Cache', 'Proxy-Authenticate')
BODY_SNIPPET_LENGTH = 500

RETRYABLE_REQUEST_EXCEPTIONS = (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError)


class RemoteResponseError(Exception):
    """ The remote, or something in the network path to it, answered with an error or with something other than JSON """

    def __init__(self, response: Response):
        self.status_code = response.status_code
        content_type = response.headers.get('Content-Type', '')
        headers = ", ".join(f"{h}={v}" for h in DIAGNOSTIC_HEADERS if (v := response.headers.get(h)))
        body = response.content[:BODY_SNIPPET_LENGTH].decode(errors='replace')
        super().__init__(f"{response.request.method} {response.url} returned {response.status_code} "
                         f"Content-Type={content_type!r} [{headers}]: {body}")

    @property
    def retryable(self) -> bool:
        # a 4xx is the remote rejecting our request (auth, bad params) - sending it again gets the same answer
        return not 400 <= self.status_code < 500


def check_json_response(response: Response) -> Response:
    content_type = response.headers.get('Content-Type', '')
    if not response.ok or not content_type.startswith('application/json'):
        raise RemoteResponseError(response)
    return response


def request_json_with_retries(send: Callable[[], Response], description: str) -> Response:
    """ send() is called again on a retry, so it must be safe to repeat (a GET, or an upsert) """
    for retries_left in range(RETRY_ATTEMPTS - 1, 0, -1):
        try:
            return check_json_response(send())
        except (*RETRYABLE_REQUEST_EXCEPTIONS, RemoteResponseError) as e:
            if isinstance(e, RemoteResponseError) and not e.retryable:
                raise
            logging.warning("%s: %s - %d retries left", description, e, retries_left)
            time.sleep(RETRY_DELAY_SECS)
    return check_json_response(send())
