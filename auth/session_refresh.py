from mozilla_django_oidc.middleware import SessionRefresh


class VariantGridSessionRefresh(SessionRefresh):

    def is_refreshable_url(self, request):
        if '/api/' in request.path or request.is_ajax():
            return False
        return super().is_refreshable_url(request)
