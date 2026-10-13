import json

from django.conf import settings
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import PermissionDenied
from django.http import HttpRequest, HttpResponseBadRequest
from django.http.response import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.utils import timezone
from django.views.decorators.http import require_POST

from eventlog.models import ApiRequestLog, Event
from library.django_utils import require_superuser
from library.enums.log_level import LogLevel


def eventlog_view(view_func):
    """ Probably don't use this as it can add 100ms of time to a page as it does a DB write """

    def _wrapped_view(request, *args, **kwargs):
        app_name = str(view_func.__module__).split('.', maxsplit=1)[0]
        Event.objects.create(user=request.user,
                             app_name=app_name,
                             name=f"viewed_{view_func.__name__}",
                             date=timezone.now(),
                             details=request.path,
                             severity=LogLevel.INFO)
        return view_func(request, *args, **kwargs)

    return _wrapped_view


def eventlog(request):
    return render(request, 'eventlog.html', context={})


def eventlog_detail(request: HttpRequest, pk: int):
    event = get_object_or_404(Event, pk=pk)
    if not (request.user.is_superuser or event.user == request.user):
        raise PermissionDenied()
    return render(request, 'eventlog_detail.html', context={'event': event})


@require_superuser
def api_request_log(request):
    view_names = ApiRequestLog.objects.order_by('view_name').values_list('view_name', flat=True).distinct()
    context = {
        "view_names": view_names,
        "view_name": request.GET.get("view_name"),
        "date_from": request.GET.get("date_from"),
        "date_to": request.GET.get("date_to"),
        "retention_days": settings.API_REQUEST_LOG_DAYS,
    }
    return render(request, 'api_request_log.html', context=context)


@require_superuser
def api_request_log_object_tab(request: HttpRequest, app_label: str, model_name: str, object_pk: str):
    """ The API calls about one object (ApiRequestLog.for_object), as a tab on its page """
    content_type = get_object_or_404(ContentType, app_label=app_label, model=model_name)
    context = {
        "app_label": app_label,
        "model_name": model_name,
        "object_pk": object_pk,
        "verbose_name": content_type.model_class()._meta.verbose_name,
    }
    return render(request, 'api_request_log_object_tab.html', context=context)


@require_superuser
def api_request_log_detail(request: HttpRequest, pk: int):
    api_request = get_object_or_404(ApiRequestLog, pk=pk)
    request_body = api_request.request_body
    try:
        request_body = json.dumps(json.loads(request_body), indent=2)
    except (TypeError, ValueError):
        pass  # not JSON, or trimmed part way through
    context = {'api_request': api_request, 'request_body': request_body}
    return render(request, 'api_request_log_detail.html', context=context)


@require_POST
def create_event(request):
    app_name = request.POST.get("app_name")
    event_name = request.POST.get("event_name")
    details = request.POST.get("details")
    severity = request.POST.get("severity")

    if not app_name or not event_name or not severity:
        return HttpResponseBadRequest("Missing required parameters")

    valid_severities = {choice[0] for choice in LogLevel.CHOICES}
    if severity not in valid_severities:
        return HttpResponseBadRequest("Invalid severity")

    event_kwargs = {
        "user": request.user,
        "app_name": app_name,
        "name": event_name,
        "date": timezone.now(),
        "details": details,
        "severity": severity,
    }
    max_lengths = {"app_name": 100, "name": 200, "details": 10_000}
    for field, max_length in max_lengths.items():
        value = event_kwargs[field]
        if value and len(value) > max_length:
            event_kwargs[field] = value[:max_length] + f"\n[trimmed to {max_length} characters]"

    Event.objects.create(**event_kwargs)
    return HttpResponse()
