import logging
import time

from django.shortcuts import redirect
from django.utils.deprecation import MiddlewareMixin
from django.utils.http import url_has_allowed_host_and_scheme
from requests import HTTPError

logger = logging.getLogger(__name__)

# A second missing-state failure within this window is left as a 400, so a broken setup can't redirect forever
MISSING_STATE_RETRY_SECONDS = 60
MISSING_STATE_RETRY_SESSION_KEY = "oidc_missing_state_retry_at"


class HandleOIDC400Middleware(MiddlewareMixin):
    """
    Turns two failures of a stale or replayed /oidc/ callback into a redirect, rather than a bare 400:

    * The callback's state isn't in the session - concurrent login flows (silent refreshes, parallel tabs) overwrite
      each other's session["oidc_states"], or the browser replays a cached callback URL (back button, reload,
      restored tabs) (mozilla-django-oidc #435, our #2348). Caught before the view raises SuspiciousOperation and
      redirected once to the page the user was going to, which starts a fresh (silent, while the Keycloak session
      lives) login.
    * Keycloak rejects a cached callback's already-used authorization code (HTTPError 400): redirect to root
      so the user goes through a fresh login.

    Any other exception is left to Django's handling
    """
    def process_view(self, request, view_func, view_args, view_kwargs):
        if self._is_unknown_callback_state(request):
            return self._retry_login(request)
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

    @staticmethod
    def _retry_login(request):
        now = time.time()
        last_retry = request.session.get(MISSING_STATE_RETRY_SESSION_KEY, 0)
        if now - last_retry < MISSING_STATE_RETRY_SECONDS:
            logger.warning("OIDC callback state missing again within %ds, not retrying", MISSING_STATE_RETRY_SECONDS)
            return None  # the view raises SuspiciousOperation

        request.session[MISSING_STATE_RETRY_SESSION_KEY] = now
        next_url = request.session.get("oidc_login_next")
        if not (next_url and url_has_allowed_host_and_scheme(next_url, allowed_hosts={request.get_host()},
                                                             require_https=request.is_secure())):
            next_url = '/'
        logger.warning("OIDC callback state not in session, retrying login via %s", next_url)
        return redirect(next_url)
