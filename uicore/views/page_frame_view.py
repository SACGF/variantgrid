"""
The page frame endpoint (#2007): JSON of uicore/page_frame.py:PageFrame for the page whose url name is passed, requested
by the first script tag of uicore/page/base.html and applied by global.js:loadPageFrame.
"""
from dataclasses import asdict

from django.contrib.auth.decorators import login_not_required
from django.http import HttpRequest, JsonResponse
from django.middleware.csrf import get_token
from django.views.decorators.cache import never_cache
from django.views.decorators.http import require_GET

from uicore.page_frame import get_page_frame


@login_not_required
@never_cache
@require_GET
def page_frame(request: HttpRequest) -> JsonResponse:
    """ Public so an expired session gets an empty page frame rather than a redirect to the login page. Sets the CSRF
        cookie, which the menu's POST links (postUrl in global.js) read now that the menu carries no token """
    get_token(request)
    return JsonResponse(asdict(get_page_frame(request, request.GET.get('url_name'))))
