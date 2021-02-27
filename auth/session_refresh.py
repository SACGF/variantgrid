from mozilla_django_oidc.middleware import SessionRefresh


class VariantGridSessionRefresh(SessionRefresh):

    def is_refreshable_url(self, request):
        if '/api/' in request.path or request.is_ajax():
            print(f"{request.path} That looks like an API, aint going to refresh check")
            return False
        print(f"{request.path} does not look like API")
        return super().is_refreshable_url(request)
