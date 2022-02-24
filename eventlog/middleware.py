from django.conf import settings
from django.urls import resolve

from eventlog.models import ViewEvent

INGORE_SEGMENTS = {"api", "datatable"}

class PageViewsMiddleware:

    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.

        response = self.get_response(request)

        # Code to be executed for each request/response after
        # the view is called.

        return response

    def process_view(self, request, view_func, view_args, view_kwargs):
        parts = request.path_info.split('/')
        if len(parts) > 1:
            if not any(segment in INGORE_SEGMENTS for segment in parts):
                app = request.path_info.split('/')[1]  # FYI not guarenteed to be the app, but closest thing I can find that links it
                if app in settings.LOG_ACTIVITY_APPS:
                    if url_obj := resolve(request.path_info):
                        # but url_obj.app_name returns an empty string

                        all_params = {**view_kwargs, **request.GET.dict(), **request.POST.dict()}

                        ViewEvent(
                            user=request.user,
                            view_name=f"{app}:{url_obj.view_name}",
                            args=all_params,
                            path=request.get_full_path(),
                            method=request.method
                        ).save()

        pass