from typing import Optional

from django.conf import settings
from oauthlib.oauth2 import LegacyApplicationClient
from requests.auth import HTTPBasicAuth
from requests_oauthlib.oauth2_auth import OAuth2
from requests_oauthlib.oauth2_session import OAuth2Session


class OAuthConnector:

    @staticmethod
    def shariant_oauth_connector(sync_details):
        if not sync_details in settings.SYNC_DETAILS:
            actual_keys = ", ".join(settings.SYNC_DETAILS.keys()) if settings.SYNC_DETAILS else "No options configured"
            raise ValueError(f"Could not find sync details '{sync_details}' options are: {actual_keys}")

        sync_details = settings.SYNC_DETAILS[sync_details]
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
                 **kwargs):
        self.username = username
        self.password = password
        self.host = host
        self.oauth_url = oauth_url
        self.client_id = client_id

        self._auth = None

    def auth(self):
        # TODO, this should be @lazy (rather than a manual implementation of @lazy)
        if not self._auth:
            if self.oauth_url:
                oauth = OAuth2Session(client=LegacyApplicationClient(client_id=self.client_id))
                token = oauth.fetch_token(token_url=self.oauth_url, username=self.username, password=self.password, client_id=self.client_id)
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
