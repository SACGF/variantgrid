import sys

from django.conf import settings
from django.core.checks import Error, Tags, register


@register(Tags.security)
def check_secure_cookies_not_used_with_runserver(app_configs, **kwargs):
    """runserver only speaks HTTP, and browsers never store a Secure cookie from an HTTP response, so
    logging in silently fails - you get sent back to the login page with no error."""
    errors = []
    if "runserver" in sys.argv:
        secure_cookie_settings = [name for name in ("SESSION_COOKIE_SECURE", "CSRF_COOKIE_SECURE")
                                  if getattr(settings, name)]
        if secure_cookie_settings:
            errors.append(
                Error(
                    f"{', '.join(secure_cookie_settings)} set while running runserver - nobody will be able to log in",
                    hint="Development machines copy variantgrid/settings/env_developers/_settings_template.py, "
                         "which leaves out the https_settings import. If you are serving HTTPS locally, silence "
                         "this with SILENCED_SYSTEM_CHECKS.",
                    id="variantgrid.E001",
                )
            )
    return errors
