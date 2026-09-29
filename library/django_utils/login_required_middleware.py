"""
Site-wide login: Django's LoginRequiredMiddleware (views opt out with django.contrib.auth.decorators.login_not_required)
plus settings.PUBLIC_PATHS, the URL prefixes that skip it. DRF authenticates at view time, so at middleware time a
token-authenticated API request is still anonymous - its prefix has to be public for DRF to answer it (or 401).
"""
import re

from django.conf import settings
from django.contrib.admin import AdminSite
from django.contrib.auth.middleware import LoginRequiredMiddleware


class PublicPathsLoginRequiredMiddleware(LoginRequiredMiddleware):
    """ Anonymous requests always go to settings.LOGIN_URL (OIDC on Shariant). Django would follow a view's own
        login_url instead, and the admin, staff_member_required and user_passes_test views name the admin's password
        login form """

    def __init__(self, get_response):
        super().__init__(get_response)
        self.public_patterns = [re.compile(public_path) for public_path in settings.PUBLIC_PATHS]

    def process_view(self, request, view_func, view_args, view_kwargs):
        if any(pattern.match(request.path) for pattern in self.public_patterns):
            return None
        # Django marks the admin login login_not_required - keep it behind the site login
        if isinstance(getattr(view_func, "__self__", None), AdminSite) and not request.user.is_authenticated:
            return self.handle_no_permission(request, view_func)
        return super().process_view(request, view_func, view_args, view_kwargs)

    def get_login_url(self, view_func):
        return settings.LOGIN_URL

    def get_redirect_field_name(self, view_func):
        return self.redirect_field_name
