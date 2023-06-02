from django import template
from django.conf import settings

register = template.Library()


@register.inclusion_tag("uicore/login/login.html", takes_context=True)
def login_form(context, next, form):

    """
    Provides a login form, or link to a 3rd party OIDC authentication provider based on settings.
    """

    return {
        "user": context.request.user,
        "next": next,
        "form": form,
        "use_oidc": settings.USE_OIDC,
        "maintenance_mode": settings.MAINTENANCE_MODE
    }
