from functools import cached_property
from typing import Optional

import requests
from django.conf import settings
from oauthlib.oauth2 import LegacyApplicationClient
from requests import Response
from requests.auth import HTTPBasicAuth
from requests_oauthlib.oauth2_auth import OAuth2
from requests_oauthlib.oauth2_session import OAuth2Session


class ServerAuth:
    """
    Previously called OAuth but handles a URL to a server with:
    BASIC_AUTH
    OAUTH
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

            oauth = OAuth2Session(client=LegacyApplicationClient(client_id=self.client_id))
            token = oauth.fetch_token(token_url=self.oauth_url, username=self.username, password=self.password, auth=auth_auth, client_id=self.client_id)
            return OAuth2(client_id=self.client_id, token=token)
        else:
            # fall back to basic auth
            return HTTPBasicAuth(
                username=self.username,
                password=self.password
            )

    def get(self, url_suffix: str, **kwargs) -> Response:
        return requests.get(
            url=self.url(url_suffix),
            auth=self.auth,
            **kwargs
        )

    def post(self, url_suffix: str, **kwargs) -> Response:
        return requests.post(
            url=self.url(url_suffix),
            auth=self.auth,
            **kwargs
        )

    def url(self, path: str):
        if path[0:1] == '/':
            path = path[1:]
        return self.host + '/' + path
