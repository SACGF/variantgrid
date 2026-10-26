import re

from mozilla_django_oidc.middleware import SessionRefresh

# DRF lives under /api/ or /<app>/api/
_API_PATH_RE = re.compile(r'^/([^/]+/)?api/')


class VariantGridSessionRefresh(SessionRefresh):

    @staticmethod
    def is_ajax(request):
        # request.is_ajax has been deprecated, is somewhat jQuery dependant
        # we need to start making sure all our ajax targets are in some kind of pattern
        return request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest'

    @staticmethod
    def is_sub_resource(request):
        """ fetch() (no X-Requested-With), iframes etc. Browsers send Sec-Fetch-Dest "document" only for top-level
            page loads; clients that don't send it are treated as page loads """
        return request.headers.get('Sec-Fetch-Dest', 'document') != 'document'

    def is_refreshable_url(self, request):
        if _API_PATH_RE.match(request.path) or self.is_ajax(request) or self.is_sub_resource(request):
            # API requests use token-based auth (OIDCAuthentication in DRF) and do not need
            # browser-based session refresh. Requests made from within a page are excluded so they don't
            # get an OIDC redirect, and so each one doesn't start a login flow writing its own state
            # to the session, racing the page load's flow.
            return False
        return super().is_refreshable_url(request)
