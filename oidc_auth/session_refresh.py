from mozilla_django_oidc.middleware import SessionRefresh


class VariantGridSessionRefresh(SessionRefresh):

    @staticmethod
    def is_ajax(request):
        # request.is_ajax has been deprecated, is somewhat jQuery dependant
        # we need to start making sure all our ajax targets are in some kind of pattern
        return request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest'

    def is_refreshable_url(self, request):
        if '/api/' in request.path or VariantGridSessionRefresh.is_ajax(request):
            return False
        return super().is_refreshable_url(request)
