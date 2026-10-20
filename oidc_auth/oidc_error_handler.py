import logging
import time

from django.core.exceptions import SuspiciousOperation
from django.shortcuts import redirect
from django.utils.deprecation import MiddlewareMixin
from django.utils.http import url_has_allowed_host_and_scheme
from requests import HTTPError

logger = logging.getLogger(__name__)

# A second missing-state failure within this window is left as a 400, so a broken setup can't redirect forever
MISSING_STATE_RETRY_SECONDS = 60
MISSING_STATE_RETRY_SESSION_KEY = "oidc_missing_state_retry_at"
# Raised by mozilla_django_oidc.views.OIDCAuthenticationCallbackView.get
MISSING_STATE_MESSAGE = "OIDC callback state not found in session"


class HandleOIDC400Middleware(MiddlewareMixin):
    """
    Turns two failures of a stale or replayed /oidc/ callback into a redirect, rather than a bare 400:

    * Keycloak rejects a cached callback's already-used authorization code (HTTPError 400): redirect to root
      so the user goes through a fresh login.
    * The callback's state isn't in the session (SuspiciousOperation) - concurrent silent refreshes overwrite
      each other's session["oidc_states"], or the callback URL is replayed (mozilla-django-oidc #435, our #2348).
      Redirect once to the page the user was going to, which starts a fresh (silent, while the Keycloak
      session lives) login.

    Any other exception is left to Django's handling
    """
    def process_exception(self, request, exception):
        if request.path.startswith('/oidc/') and self._is_token_rejected(exception):
            logger.warning("OIDC provider rejected token request, redirecting to root: %s", exception)
            return redirect('/')
        if self._is_missing_state(request, exception):
            return self._retry_login(request)
        return None

    @staticmethod
    def _is_token_rejected(exception) -> bool:
        return (isinstance(exception, HTTPError) and exception.response is not None
                and exception.response.status_code == 400)

    @staticmethod
    def _is_missing_state(request, exception) -> bool:
        # The callback also raises SuspiciousOperation from authenticate() (bad claims, nonce) after it has
        # consumed the state, so match the library's message rather than checking the session
        return (request.path == '/oidc/callback/' and isinstance(exception, SuspiciousOperation)
                and str(exception).startswith(MISSING_STATE_MESSAGE))

    @staticmethod
    def _retry_login(request):
        now = time.time()
        last_retry = request.session.get(MISSING_STATE_RETRY_SESSION_KEY, 0)
        if now - last_retry < MISSING_STATE_RETRY_SECONDS:
            logger.warning("OIDC callback state missing again within %ds, not retrying", MISSING_STATE_RETRY_SECONDS)
            return None

        request.session[MISSING_STATE_RETRY_SESSION_KEY] = now
        next_url = request.session.get("oidc_login_next")
        if not (next_url and url_has_allowed_host_and_scheme(next_url, allowed_hosts={request.get_host()},
                                                             require_https=request.is_secure())):
            next_url = '/'
        logger.warning("OIDC callback state not in session, retrying login via %s", next_url)
        return redirect(next_url)
