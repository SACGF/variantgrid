import json
import time
from typing import Any, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.http.response import HttpResponseRedirectBase
from django.urls import Resolver404, resolve

from eventlog.models import (
    API_KEY_PREFIX,
    ApiRequestLog,
    IntegrationActivity,
    IntegrationActivityStatus,
    ViewEvent,
)
from library.integration_status import IntegrationDirection
from oidc_auth.session_refresh import VariantGridSessionRefresh

IGNORE_SEGMENTS = {"api", "datatable", "citations_json"}  # this should be mostly redundant to is_ajax call
IGNORE_TEXT = {"detail", "metrics"}
IGNORE_VIEW_NAME_SUFFIX = {"_detail", "_autocomplete"}
IGNORE_PARAMETERS = {"csrfmiddlewaretoken"}

REDACTED = "[redacted]"
REDACTED_KEY_PARTS = ("password", "token", "secret", "api_key")


class PageViewsMiddleware:

    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.

        response = self.get_response(request)
        # don't record redirects, as the view that we're redirected to will be enough
        if hasattr(request, 'view_event'):
            view_event: Optional[ViewEvent]
            if view_event := request.view_event:
                # if a user got taken directly to something from a search, record it
                if not isinstance(response, HttpResponseRedirectBase) or view_event.view_name == 'variantopedia:search':
                    view_event.save()

        # Code to be executed for each request/response after
        # the view is called.

        return response

    # This code adds too much complexity and redundancy to the logging ()
    # leave classification activity to the

    # def is_special_case(self, request, view_func, view_args, view_kwargs) -> Optional[ViewEvent]:
    #     if 'classification' not in settings.LOG_ACTIVITY_APPS:
    #         return
    #
    #     # check to HTTP_X_REQUESTED_WITH as that implies it's a browser performing the request so hopefully we skip it when it's being performed by the API
    #     if request.path_info == '/classification/api/classifications/v1/record/' and request.method == 'POST' and request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest':
    #         data = json.loads(request.body)
    #         if isinstance(data, dict) and data.get('source') == 'form':
    #             all_params = {'classification_id': data.get('id').get('record_id'), **data.get('patch')}
    #             if data.get('delete') is not None:
    #                 all_params['delete'] = data.get('delete')
    #             if data.get('publish') is not None:
    #                 all_params['publish'] = data.get('publish')
    #
    #             return ViewEvent(
    #                 user=request.user,
    #                 # doesn't match any actual view name, just trying to be explicit
    #                 view_name='classification:classification_form',
    #                 args=all_params,
    #                 path=request.get_full_path(),
    #                 method=request.method,
    #                 referer=request.headers.get('Referer')
    #             )
    #
    #         pass

    def process_view(self, request, view_func, view_args, view_kwargs):
        # if special_view_event := self.is_special_case(request, view_func, view_args, view_kwargs):
        #     request.view_event = special_view_event
        #     return

        if VariantGridSessionRefresh.is_ajax(request):
            return

        parts = request.path_info.split('/')
        if len(parts) > 1:
            if any(ignore_text in request.path_info for ignore_text in IGNORE_TEXT):
                return

            if not any(segment in IGNORE_SEGMENTS for segment in parts):
                app = request.path_info.split('/')[1]  # FYI not guaranteed to be the app, but closest thing I can find that links it
                if app in settings.LOG_ACTIVITY_APPS:
                    if url_obj := resolve(request.path_info):
                        view_name = url_obj.view_name
                        if any(view_name.endswith(ignore_suffix) for ignore_suffix in IGNORE_VIEW_NAME_SUFFIX):
                            return

                        # but url_obj.app_name returns an empty string

                        all_params = {**view_kwargs}
                        for key, value in {**request.GET.dict(), **request.POST.dict()}.items():
                            key: str
                            value: Any
                            if value.lower() == "true":
                                value = True
                            elif value.lower() == "false":
                                value = False
                            else:
                                try:
                                    value = int(value)
                                except ValueError:
                                    pass
                            all_params[key] = value

                        for ignore_me in IGNORE_PARAMETERS:
                            if ignore_me in all_params:
                                all_params.pop(ignore_me)

                        if app == "ontology" and (term := all_params.get("term")):
                            # ontology has to replace : with _ so as to not confuse URLs
                            all_params["term"] = term.replace("_", ":")

                        # hack to split up classification ID when it's in the form of "classification_id.modification_timestamp"
                        if classification_id := all_params.get("classification_id"):
                            if isinstance(classification_id, (str, float)):
                                check_for_parts = str(classification_id)
                                parts = check_for_parts.split(".")
                                try:
                                    all_params["classification_id"] = int(parts[0])
                                    if len(parts) > 1:
                                        all_params["modification_timestamp"] = float(parts[1])
                                except ValueError:
                                    pass

                        # call then later captures this and saves it to the DB (if appropriate)
                        user: Optional[User] = None
                        if request_user := request.user:
                            if request_user.is_authenticated:
                                user = request_user

                        if user:
                            request.view_event = ViewEvent(
                                user=user,
                                view_name=f"{app}:{url_obj.view_name}",
                                args=all_params,
                                path=request.get_full_path(),
                                method=request.method,
                                referer=request.headers.get('Referer')
                            )


class IntegrationApiMiddleware:
    """
    Records each call an external client makes to us: the rolling "last run" on Server Status
    (IntegrationActivity) and one ApiRequestLog row with what was sent and what came back, for
    endpoints that would otherwise leave nothing behind but the records they wrote.

    Driven by settings.INTEGRATION_API_TRACKING - {url path prefix: display name} - so a
    deployment only pays for the endpoints it names.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        tracking = settings.INTEGRATION_API_TRACKING
        display_name = next((name for prefix, name in tracking.items() if request.path_info.startswith(prefix)), None)
        if display_name is None:
            return self.get_response(request)

        # read before the view, which consumes the stream
        request_body = _request_body_for_log(request)
        start = time.perf_counter()
        response = self.get_response(request)
        duration_ms = int((time.perf_counter() - start) * 1000)
        self._record(request, response, display_name, request_body, duration_ms)
        return response

    @staticmethod
    def _record(request, response, display_name: str, request_body: Optional[str], duration_ms: int):
        try:
            resolver_match = request.resolver_match or resolve(request.path_info)
        except Resolver404:
            return
        view_name = resolver_match.view_name

        status_code = response.status_code
        if status_code >= 400:
            status = IntegrationActivityStatus.ERROR
            message = f"HTTP {status_code} from {request.path_info}"
        else:
            status = IntegrationActivityStatus.SUCCESS
            message = None

        IntegrationActivity.record(
            key=f"{API_KEY_PREFIX}{view_name}",
            name=f"{display_name} ({view_name})",
            direction=IntegrationDirection.INBOUND,
            status=status,
            changed=status == IntegrationActivityStatus.SUCCESS and request.method != "GET",
            message=message
        )

        response_body = None
        if status_code >= 400 and not response.streaming:
            response_body = _trim_for_log(response.content.decode(errors="replace"))

        # A view names the model its calls are about with api_log_model, and the pks it found with
        # log_api_request_objects - a call that fails before finding one keeps only the model
        view_class = getattr(resolver_match.func, "cls", None) or getattr(resolver_match.func, "view_class", None)
        content_type = None
        if api_log_model := getattr(view_class, "api_log_model", None):
            content_type = ContentType.objects.get_for_model(api_log_model)

        user = getattr(request, "user", None)  # DRF copies its authenticated user back onto the Django request
        ApiRequestLog.objects.create(
            user=user if user and user.is_authenticated else None,
            method=request.method,
            path=request.path_info,
            query_string=request.META.get("QUERY_STRING", ""),
            view_name=view_name,
            view_kwargs=resolver_match.kwargs,
            status_code=status_code,
            duration_ms=duration_ms,
            request_body=request_body,
            response_body=response_body,
            content_type=content_type,
            object_pks=sorted(getattr(request, "api_log_object_pks", ())),
        )


def _trim_for_log(text: str) -> str:
    max_length = settings.API_REQUEST_LOG_MAX_BODY
    if len(text) > max_length:
        text = text[:max_length] + f"\n[trimmed from {len(text)} characters]"
    return text


def _redact(data):
    if isinstance(data, dict):
        return {key: REDACTED if any(part in str(key).lower() for part in REDACTED_KEY_PARTS) else _redact(value)
                for key, value in data.items()}
    if isinstance(data, list):
        return [_redact(value) for value in data]
    return data


def _request_body_for_log(request) -> Optional[str]:
    """ Uploads and big bodies are noted by size rather than read, so they stay out of memory and
        under DATA_UPLOAD_MAX_MEMORY_SIZE, which request.body enforces and DRF's stream parsing doesn't """
    try:
        content_length = int(request.META.get("CONTENT_LENGTH") or 0)
    except ValueError:
        content_length = 0
    if not content_length:
        return None

    content_type = request.content_type or "unknown content type"
    max_read = min(filter(None, [settings.API_REQUEST_LOG_MAX_READ_BYTES, settings.DATA_UPLOAD_MAX_MEMORY_SIZE]))
    if content_type.startswith("multipart/") or content_length > max_read:
        return f"[{content_type}, {content_length} bytes]"

    body = request.body.decode(errors="replace")
    if content_type == "application/json":
        try:
            body = json.dumps(_redact(json.loads(body)))
        except ValueError:
            pass
    return _trim_for_log(body)
