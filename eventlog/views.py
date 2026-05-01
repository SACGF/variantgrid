from django.core.exceptions import PermissionDenied
from django.http import HttpRequest, HttpResponseBadRequest
from django.http.response import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.utils import timezone
from django.views.decorators.http import require_POST

from eventlog.models import Event
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

    app_name = app_name[:100]
    event_name = event_name[:200]
    if details:
        details = details[:10_000]

    Event.objects.create(user=request.user,
                         app_name=app_name,
                         name=event_name,
                         date=timezone.now(),
                         details=details,
                         severity=severity)
    return HttpResponse()
