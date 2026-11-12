"""
API token pages for snpdb/models/models_api_token.py:ApiToken: a user's own tokens (api_tokens - create, revoke, use
history) and the admin list of every user's tokens (api_tokens_admin). The user pages follow
settings.API_TOKEN_SELF_SERVICE_ENABLED through URLS_NAME_REGISTER_SETTINGS, so they are admin-only where it is off.
"""
from datetime import timedelta
from urllib.parse import urlencode

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import Count, Q, QuerySet
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils import timezone
from django.views.decorators.http import require_POST

from library.django_utils import require_superuser
from snpdb.forms import ApiTokenCreateForm
from snpdb.models.models_api_token import ApiToken

RECENT_USE_DAYS = 30
MAX_USES_SHOWN = 500


def _with_recent_uses(qs: QuerySet[ApiToken]) -> QuerySet[ApiToken]:
    since = timezone.now() - timedelta(days=RECENT_USE_DAYS)
    return qs.select_related("user", "revoked_by").annotate(
        recent_uses=Count("uses", filter=Q(uses__used__gte=since))).order_by("-created")


def _get_token(user: User, digest: str) -> ApiToken:
    token = get_object_or_404(ApiToken.objects.select_related("user"), digest=digest)
    if not token.can_manage(user):
        raise PermissionDenied("Not your API token")
    return token


def api_tokens(request):
    form = ApiTokenCreateForm(request.POST or None)
    new_token = new_token_key = None
    if form.is_valid():
        new_token, new_token_key = ApiToken.create_for_user(request.user, form.cleaned_data["name"], form.expiry,
                                                            created_by=request.user)
        form = ApiTokenCreateForm()

    context = {
        "form": form,
        "new_token": new_token,
        "new_token_key": new_token_key,
        "tokens": _with_recent_uses(ApiToken.objects.filter(user=request.user)),
        "recent_use_days": RECENT_USE_DAYS,
    }
    return render(request, "snpdb/settings/api_tokens.html", context)


@require_superuser
def api_tokens_admin(request):
    context = {
        "tokens": _with_recent_uses(ApiToken.objects.all()),
        "recent_use_days": RECENT_USE_DAYS,
    }
    return render(request, "snpdb/settings/api_tokens_admin.html", context)


def api_token_detail(request, digest: str):
    token = _get_token(request.user, digest)
    api_log_url = None
    if request.user.is_superuser and settings.INTEGRATION_API_TRACKING:
        api_log_url = reverse("api_request_log") + "?" + urlencode({"auth_token": token.pk})
    context = {
        "token": token,
        "uses": token.uses.order_by("-used")[:MAX_USES_SHOWN],
        "max_uses_shown": MAX_USES_SHOWN,
        "api_log_url": api_log_url,
    }
    return render(request, "snpdb/settings/api_token_detail.html", context)


@require_POST
def api_token_revoke(request):
    token = _get_token(request.user, request.POST["digest"])
    token.revoke(request.user)
    messages.add_message(request, messages.INFO, f"Revoked API token '{token.name}'")
    if request.POST.get("admin") and request.user.is_superuser:
        return redirect("api_tokens_admin")
    return redirect("api_tokens")
