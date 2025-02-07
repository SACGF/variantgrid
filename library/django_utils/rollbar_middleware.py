from rollbar.contrib.django.middleware import RollbarNotifierMiddleware


class RollbarIgnoreException(Exception):
    """ Throw an error that subclasses this for Rollbar to ignore """


class CustomRollbarNotifierMiddleware(RollbarNotifierMiddleware):

    def get_extra_data(self, request, exc):
        if isinstance(exc, Exception):
            return {"exception_message": str(exc)}
        return None

    def get_payload_data(self, request, exc):
        return

    def process_exception(self, request, exc):
        """ Copied from rollbar.contrib.django.middleware.RollbarNotifierMiddlewareExcluding404 """
        if isinstance(exc, RollbarIgnoreException):
            return
        else:
            super().process_exception(request, exc)
