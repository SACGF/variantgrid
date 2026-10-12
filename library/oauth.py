import os
from functools import cached_property
from typing import Optional
from urllib.parse import urljoin, urlparse

import requests
from django.conf import settings
from oauthlib.oauth2 import LegacyApplicationClient
from requests import Response
from requests.auth import HTTPBasicAuth
from requests_oauthlib.oauth2_auth import OAuth2
from requests_oauthlib.oauth2_session import OAuth2Session

from library.constants import MINUTE_SECS

# Keycloak grants its default client scopes (e.g. "email profile") on top of the requested
# "openid" - stop oauthlib raising "Scope has changed" over the widened grant
os.environ.setdefault("OAUTHLIB_RELAX_TOKEN_SCOPE", "1")


class ServerAuth:
    """
    Previously called OAuth but handles a URL to a server with:
    BASIC_AUTH
    OAUTH and
    OAUTH where you need BASIC_AUTH first
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
                 **kwargs):
        self.username = username
        self.password = password
        self.host = host
        self.oauth_url = oauth_url
        self.client_id = client_id
        self.app_username = app_username
        self.app_password = app_password

    @cached_property
    def auth(self):
        if self.oauth_url:
            auth_auth = None
            if self.app_username or self.app_password:
                auth_auth = HTTPBasicAuth(
                    username=self.app_username,
                    password=self.app_password
                )

            # scope must be set on the oauthlib client (not OAuth2Session) for it to be
            # included in the token request body - Keycloak requires "openid"
            oauth = OAuth2Session(client=LegacyApplicationClient(client_id=self.client_id, scope="openid"))
            # include_client_id sends client_id in the request body (public client) - without
            # it, requests_oauthlib fabricates a Basic auth header with an empty client_secret,
            # which Keycloak rejects with invalid_client
            token = oauth.fetch_token(token_url=self.oauth_url, username=self.username, password=self.password, auth=auth_auth, include_client_id=True)
            return OAuth2(client_id=self.client_id, token=token)
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
        scheme = urlparse(self.host).scheme
        if scheme not in ('https', 'http'):
            raise ValueError(f"ServerAuth host must use http(s), got scheme {scheme!r}")
        base = self.host if self.host.endswith('/') else self.host + '/'
        return urljoin(base, path.lstrip('/'))
