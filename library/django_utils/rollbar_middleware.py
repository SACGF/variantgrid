from rollbar.contrib.django.middleware import RollbarNotifierMiddleware


class RollbarIgnoreException(Exception):
    """ Throw an error that subclasses this for Rollbar to ignore """


class CustomRollbarNotifierMiddleware(RollbarNotifierMiddleware):

    def get_extra_data(self, request, exc):
        if isinstance(exc, Exception):
            return {"exception_message": str(exc)}
        return

    def get_payload_data(self, request, exc):
        return

    """ Copied from rollbar.contrib.django.middleware.RollbarNotifierMiddlewareExcluding404 """
    def process_exception(self, request, exc):
        if isinstance(exc, RollbarIgnoreException):
            return
        else:
            super().process_exception(request, exc)
