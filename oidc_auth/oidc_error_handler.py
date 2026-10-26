import logging

from django.shortcuts import redirect
from django.utils.deprecation import MiddlewareMixin
from requests import HTTPError

logger = logging.getLogger(__name__)


class HandleOIDC400Middleware(MiddlewareMixin):
    """
    Browsers often replay a cached /oidc/ callback URL (back button, reload, restored tabs), and parallel login flows
    overwrite each other's state in the session. Both are redirected to root so the user goes through a fresh login
    (or straight to the page if already logged in), instead of a 400 page:

    * a callback whose state is no longer in the session - mozilla_django_oidc would raise SuspiciousOperation
    * Keycloak rejecting a stale authorization code with a 400

    Any other exception is left to Django's handling
    """
    def process_view(self, request, view_func, view_args, view_kwargs):
        if self._is_unknown_callback_state(request):
            logger.warning("OIDC callback state not in session (replayed callback or parallel login), redirecting to root")
            return redirect('/')
        return None

    def process_exception(self, request, exception):
        if request.path.startswith('/oidc/') and self._is_token_rejected(exception):
            logger.warning("OIDC provider rejected token request, redirecting to root: %s", exception)
            return redirect('/')
        return None

    @staticmethod
    def _is_unknown_callback_state(request) -> bool:
        """ The condition under which OIDCAuthenticationCallbackView raises SuspiciousOperation.
            No "oidc_states" at all is left to the view, which treats it as a login failure """
        if request.resolver_match.url_name != "oidc_authentication_callback":
            return False
        state = request.GET.get("state")
        oidc_states = request.session.get("oidc_states")
        return bool(state and "code" in request.GET and oidc_states is not None and state not in oidc_states)

    @staticmethod
    def _is_token_rejected(exception) -> bool:
        return (isinstance(exception, HTTPError) and exception.response is not None
                and exception.response.status_code == 400)
