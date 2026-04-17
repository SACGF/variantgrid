import logging

from django.shortcuts import redirect
from django.utils.deprecation import MiddlewareMixin

logger = logging.getLogger(__name__)


class HandleOIDC400Middleware(MiddlewareMixin):
    """
    Very frequently we're getting exceptions from talking to Rollbar with an old token (looks like it's just browser
    caching that URL - so if we get a 400 from /oidc/ redirect to root and things solve themselves - as we just send
    up to date token)
    """
    def process_exception(self, request, exception):
        # if (isinstance(exception, HTTPError) and exception.response.status_code == 400)
        # or isinstance(exception, SuspiciousActivity):
        if request.path.startswith('/oidc/'):
            logger.exception("OIDC exception at %s", request.path)
            return redirect('/')
        return None
