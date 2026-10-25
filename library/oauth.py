import os
import time
from functools import cached_property
from typing import Optional
from urllib.parse import urljoin, urlparse

import requests
from django.conf import settings
from oauthlib.oauth2 import LegacyApplicationClient, OAuth2Error
from requests import Response
from requests.auth import HTTPBasicAuth
from requests_oauthlib.oauth2_auth import OAuth2
from requests_oauthlib.oauth2_session import OAuth2Session

from library.constants import MINUTE_SECS

# Keycloak grants its default client scopes (e.g. "email profile") on top of the requested
# "openid" - stop oauthlib raising "Scope has changed" over the widened grant
os.environ.setdefault("OAUTHLIB_RELAX_TOKEN_SCOPE", "1")

LOOPBACK_HOSTNAMES = {"localhost", "127.0.0.1", "::1"}

# refresh this long before the access token expires, so it doesn't lapse between being sent and being checked
TOKEN_EXPIRY_MARGIN_SECS = 30


def url_scheme_problem(url: str, allow_http: bool = False) -> Optional[str]:
    """ Why credentials can't be sent to url, or None: only https, or http to this machine unless allow_http """
    parsed = urlparse(url)
    if parsed.scheme == 'https':
        return None
    if parsed.scheme == 'http':
        if allow_http or parsed.hostname in LOOPBACK_HOSTNAMES:
            return None
        return f"{url!r} uses http - use https, or set allow_http in the sync details"
    return f"{url!r} must use http(s), got scheme {parsed.scheme!r}"


class ServerAuth:
    """
    Previously called OAuth but handles a URL to a server with:
    BASIC_AUTH
    OAUTH and
    OAUTH where you need BASIC_AUTH first

    Credentials only go over https, or http to this machine - other http hosts need allow_http in the sync details

    OAUTH logs in once per ServerAuth, then refreshes the token as it expires (a new login if the refresh is refused),
    so share one instance across a run's requests
    """

    @staticmethod
    def for_sync_details(sync_details):
        return ServerAuth(**sync_details)

    @staticmethod
    def keycloak_connector():
        sync_details = settings.KEYCLOAK_SYNC_DETAILS
        return ServerAuth(**sync_details)

    def __init__(self,
                 username: str,
                 password: str,
                 host: str,
                 oauth_url: Optional[str] = None,
                 client_id: Optional[str] = None,
                 app_username: Optional[str] = None,
                 app_password: Optional[str] = None,
                 allow_http: bool = False,
                 **kwargs):
        self.username = username
        self.password = password
        self.host = host
        self.oauth_url = oauth_url
        self.client_id = client_id
        self.app_username = app_username
        self.app_password = app_password
        self.allow_http = allow_http

    def _check_scheme(self, url: str):
        if problem := url_scheme_problem(url, self.allow_http):
            raise ValueError(f"ServerAuth {problem}")

    @cached_property
    def _app_auth(self) -> Optional[HTTPBasicAuth]:
        if self.app_username or self.app_password:
            return HTTPBasicAuth(username=self.app_username, password=self.app_password)
        return None

    def _fetch_token(self, session: OAuth2Session):
        # include_client_id sends client_id in the request body (public client) - without
        # it, requests_oauthlib fabricates a Basic auth header with an empty client_secret,
        # which Keycloak rejects with invalid_client
        session.fetch_token(token_url=self.oauth_url, username=self.username, password=self.password,
                            auth=self._app_auth, include_client_id=True)

    @cached_property
    def _oauth_session(self) -> OAuth2Session:
        self._check_scheme(self.oauth_url)
        # scope must be set on the oauthlib client (not OAuth2Session) for it to be
        # included in the token request body - Keycloak requires "openid"
        session = OAuth2Session(client=LegacyApplicationClient(client_id=self.client_id, scope="openid"),
                                auto_refresh_kwargs={"client_id": self.client_id})
        self._fetch_token(session)
        return session

    def _refresh_if_expiring(self, session: OAuth2Session):
        # refreshed here rather than by OAuth2Session.request's auto_refresh_url, which passes the request's own
        # kwargs (json, params, stream) into the refresh request body
        expires_at = session.token.get("expires_at")
        if expires_at is None or expires_at - TOKEN_EXPIRY_MARGIN_SECS > time.time():
            return
        if session.token.get("refresh_token"):
            try:
                session.refresh_token(self.oauth_url, auth=self._app_auth)
                return
            except OAuth2Error:
                pass  # refresh token expired or revoked
        self._fetch_token(session)

    @property
    def auth(self):
        if self.oauth_url:
            session = self._oauth_session
            self._refresh_if_expiring(session)
            return OAuth2(client_id=self.client_id, token=session.token)
        else:
            # fall back to basic auth
            return HTTPBasicAuth(
                username=self.username,
                password=self.password
            )

    def get(self, url_suffix: str, timeout: int = MINUTE_SECS, **kwargs) -> Response:
        return requests.get(
            url=self.url(url_suffix),
            auth=self.auth,
            timeout=timeout,
            **kwargs
        )

    def post(self, url_suffix: str, timeout: int = MINUTE_SECS, **kwargs) -> Response:
        return requests.post(
            url=self.url(url_suffix),
            auth=self.auth,
            timeout=timeout,
            **kwargs
        )

    def url(self, path: str) -> str:
        """ path is relative to host even with a leading '/', so a host with a sub-path keeps it """
        self._check_scheme(self.host)
        base = self.host if self.host.endswith('/') else self.host + '/'
        return urljoin(base, path.lstrip('/'))
