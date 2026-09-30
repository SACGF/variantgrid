"""
The request being handled, for code that isn't passed one: logging (`library/log_utils.py`), per-request caches
(get_request_variable / set_request_variable, eg UserSettingsManager) and the current user's genome build.

RequestContextMiddleware sets the request, and clears it when the response is closed - after any streaming
content has been sent, so streamed exports still see it. Being a ContextVar rather than a thread local, each
request (or asyncio task) sees its own. Outside a request (celery, management commands, tests calling code
directly) get_current_request() is None, and request variables fall back to a store that lasts as long as the
thread's context.
"""
from contextvars import ContextVar
from typing import Any, Optional

from django.core.signals import request_finished
from django.http import HttpRequest

_current_request: ContextVar[Optional[HttpRequest]] = ContextVar("current_request", default=None)
_variables_without_request: ContextVar[Optional[dict[str, Any]]] = ContextVar("variables_without_request",
                                                                              default=None)


def get_current_request() -> Optional[HttpRequest]:
    return _current_request.get()


def get_current_user():
    if request := get_current_request():
        return getattr(request, "user", None)
    return None


def _get_variables_without_request() -> dict[str, Any]:
    if (variables := _variables_without_request.get()) is None:
        variables = {}
        _variables_without_request.set(variables)
    return variables


def set_thread_variable(key, val):
    """ 'request' sets (or with None, clears) the current request, eg to drop one a test client left behind """
    if key == "request":
        _current_request.set(val)
    else:
        _get_variables_without_request()[key] = val


def get_thread_variable(key, default=None):
    if key == "request":
        return get_current_request()
    return _get_variables_without_request().get(key, default)


def _get_request_variables(request: HttpRequest) -> dict[str, Any]:
    if (variables := getattr(request, "_variables", None)) is None:
        variables = {}
        request._variables = variables
    return variables


def get_request_variable(key, default=None, use_threadlocal_if_no_request: bool = True):
    """ With no current request, reads the store that outlives requests - or with use_threadlocal_if_no_request=False
        raises RuntimeError, for a cache that must never outlive a request (QuerySetRequestCache) """
    if request := get_current_request():
        return _get_request_variables(request).get(key, default)
    if not use_threadlocal_if_no_request:
        raise RuntimeError("No current request - is RequestContextMiddleware installed?")
    return get_thread_variable(key, default)


def set_request_variable(key, val, use_threadlocal_if_no_request: bool = True):
    if request := get_current_request():
        _get_request_variables(request)[key] = val
    elif use_threadlocal_if_no_request:
        set_thread_variable(key, val)
    else:
        raise RuntimeError("No current request - is RequestContextMiddleware installed?")


class RequestContextMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        _current_request.set(request)
        return self.get_response(request)


def _clear_current_request(**kwargs):
    _current_request.set(None)


# Sent when the response is closed, ie after streaming content has been consumed
request_finished.connect(_clear_current_request, dispatch_uid="library.request_context.clear_current_request")
