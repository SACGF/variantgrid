from typing import Optional

from django.conf import settings
from oauthlib.oauth2 import LegacyApplicationClient
from requests.auth import HTTPBasicAuth
from requests_oauthlib.oauth2_auth import OAuth2
from requests_oauthlib.oauth2_session import OAuth2Session


class OAuthConnector:

    @staticmethod
    def shariant_oauth_connector(sync_details):
        return OAuthConnector(**sync_details)

    @staticmethod
    def keycloak_connector():
        sync_details = settings.KEYCLOAK_SYNC_DETAILS
        return OAuthConnector(**sync_details)

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

        self._auth = None

    def auth(self):
        # TODO, this should be @cached_property (rather than a manual implementation of @cached_property)
        if not self._auth:
            if self.oauth_url:
                auth_auth = None
                if self.app_username or self.app_password:
                    auth_auth = HTTPBasicAuth(
                        username=self.app_username,
                        password=self.app_password
                    )

                oauth = OAuth2Session(client=LegacyApplicationClient(client_id=self.client_id))
                token = oauth.fetch_token(token_url=self.oauth_url, username=self.username, password=self.password, auth=auth_auth, client_id=self.client_id)
                self._auth = OAuth2(client_id=self.client_id, token=token)
            else:
                # fall back to basic auth
                self._auth = HTTPBasicAuth(
                    username=self.username,
                    password=self.password
                )
        return self._auth

    def url(self, path: str):
        if path[0:1] == '/':
            path = path[1:]
        return self.host + '/' + path
