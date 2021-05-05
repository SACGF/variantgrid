import sys

from django.http import Http404
from rollbar.contrib.django.middleware import RollbarNotifierMiddleware


class RollbarIgnoreException(Exception):
    """ Throw an error that subclasses this for Rollbar to ignore """
    pass


class CustomRollbarNotifierMiddleware(RollbarNotifierMiddleware):
    """ Copied from rollbar.contrib.django.middleware.RollbarNotifierMiddlewareExcluding404 """
    def process_exception(self, request, exc):
        if isinstance(exc, Http404):
            request._rollbar_notifier_original_http404_exc_info = sys.exc_info()
        elif isinstance(exc, RollbarIgnoreException):
            return
        else:
            super(self).process_exception(request, exc)