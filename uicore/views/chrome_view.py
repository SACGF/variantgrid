"""
The page chrome endpoint (#2007): JSON of uicore/chrome.py:PageChrome for the page whose url name is passed, fetched
by global.js:loadPageChrome on every page that extends uicore/page/base.html.
"""
from dataclasses import asdict

from django.contrib.auth.decorators import login_not_required
from django.http import HttpRequest, JsonResponse
from django.middleware.csrf import get_token
from django.views.decorators.cache import never_cache
from django.views.decorators.http import require_GET

from uicore.chrome import get_page_chrome


@login_not_required
@never_cache
@require_GET
def page_chrome(request: HttpRequest) -> JsonResponse:
    """ Public so an expired session gets an empty chrome rather than a redirect to the login page. Sets the CSRF
        cookie, which the menu's POST links (postUrl in global.js) read now that the menu carries no token """
    get_token(request)
    return JsonResponse(asdict(get_page_chrome(request, request.GET.get('url_name'))))
