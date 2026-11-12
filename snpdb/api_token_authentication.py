"""
DRF authentication for snpdb/models/models_api_token.py:ApiToken - knox's TokenAuthentication, listed in
REST_FRAMEWORK's DEFAULT_AUTHENTICATION_CLASSES ahead of DRF's own TokenAuthentication.

Both classes read the `Authorization: Token <key>` header and refuse a key they don't hold rather than passing it on,
so this one only takes keys with knox's TOKEN_PREFIX and leaves the rest (legacy rest_framework.authtoken keys,
40 hex characters) to DRF's. ApiTokenScheme describes it in the OpenAPI schema.
"""
from axes.helpers import get_client_ip_address
from drf_spectacular.extensions import OpenApiAuthenticationExtension
from drf_spectacular.plumbing import build_bearer_security_scheme_object
from knox.auth import TokenAuthentication
from knox.settings import knox_settings
from rest_framework.authentication import get_authorization_header

from snpdb.models.models_api_token import ApiToken


class ApiTokenAuthentication(TokenAuthentication):

    def authenticate(self, request):
        auth = get_authorization_header(request).split()
        if len(auth) == 2 and not auth[1].startswith(knox_settings.TOKEN_PREFIX.encode()):
            return None
        if result := super().authenticate(request):
            token: ApiToken = result[1]
            token.record_use(get_client_ip_address(request), request.path)
        return result

    def _cleanup_token(self, auth_token: ApiToken) -> bool:
        """ Knox deletes expired tokens here (and the user's other expired ones); ours are kept for their history,
            so this only says whether to skip the token """
        return not auth_token.is_valid


class ApiTokenScheme(OpenApiAuthenticationExtension):
    target_class = ApiTokenAuthentication
    name = 'apiTokenAuth'

    def get_security_definition(self, auto_schema):
        scheme = build_bearer_security_scheme_object(header_name='Authorization',
                                                     token_prefix=knox_settings.AUTH_HEADER_PREFIX)
        scheme["description"] = f"An API token ({knox_settings.TOKEN_PREFIX}...) from Settings > API Tokens, " \
                                f"sent as \"{knox_settings.AUTH_HEADER_PREFIX} <key>\""
        return scheme
