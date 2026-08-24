from typing import Any, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.http.response import HttpResponseRedirectBase
from django.urls import Resolver404, resolve

from eventlog.models import (
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
    Records when each external client last called us - the API equivalent of a "last run" log,
    for endpoints that would otherwise leave nothing behind but the records they wrote.

    Driven by settings.INTEGRATION_API_TRACKING - {url path prefix: display name} - so a
    deployment only pays for the endpoints it names.
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)

        tracking = settings.INTEGRATION_API_TRACKING
        if not tracking:
            return response

        for prefix, display_name in tracking.items():
            if request.path_info.startswith(prefix):
                self._record(request, response, display_name)
                break

        return response

    @staticmethod
    def _record(request, response, display_name: str):
        try:
            view_name = resolve(request.path_info).view_name
        except Resolver404:
            return

        status_code = response.status_code
        if status_code >= 400:
            status = IntegrationActivityStatus.ERROR
            message = f"HTTP {status_code} from {request.path_info}"
        else:
            status = IntegrationActivityStatus.SUCCESS
            message = None

        IntegrationActivity.record(
            key=f"api:{view_name}",
            name=f"{display_name} ({view_name})",
            direction=IntegrationDirection.INBOUND,
            status=status,
            changed=status == IntegrationActivityStatus.SUCCESS and request.method != "GET",
            message=message
        )
