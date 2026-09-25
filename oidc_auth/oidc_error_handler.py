import logging

from django.shortcuts import redirect
from django.utils.deprecation import MiddlewareMixin
from requests import HTTPError

logger = logging.getLogger(__name__)


class HandleOIDC400Middleware(MiddlewareMixin):
    """
    Browsers often replay a cached /oidc/ callback URL, so Keycloak rejects its stale authorization code with a 400.
    Redirect to root so the user goes through a fresh login. Any other exception is left to Django's handling
    """
    def process_exception(self, request, exception):
        if request.path.startswith('/oidc/') and self._is_token_rejected(exception):
            logger.warning("OIDC provider rejected token request, redirecting to root: %s", exception)
            return redirect('/')
        return None

    @staticmethod
    def _is_token_rejected(exception) -> bool:
        return (isinstance(exception, HTTPError) and exception.response is not None
                and exception.response.status_code == 400)
