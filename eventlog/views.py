from django.http.response import HttpResponse
from django.shortcuts import render
from django.utils import timezone
from django.views.decorators.http import require_POST
import logging

from eventlog.grids import EventColumns
from eventlog.models import Event
from library.enums.log_level import LogLevel


def eventlog_view(view_func):
    """ Probably don't use this as it can add 100ms of time to a page as it does a DB write """

    def _wrapped_view(request, *args, **kwargs):
        app_name = str(view_func.__module__).split('.')[0]
        Event.objects.create(user=request.user,
                             app_name=app_name,
                             name=f"viewed_{view_func.__name__}",
                             date=timezone.now(),
                             details=request.path,
                             severity=LogLevel.INFO)
        return view_func(request, *args, **kwargs)

    return _wrapped_view


def eventlog(request):
    return render(request, 'eventlog.html', context={
        'datatable_config': EventColumns(request)
    })


@require_POST
def create_event(request):
    logging.info("create_event POST=%s", str(request.POST))
    app_name = request.POST["app_name"]
    event_name = request.POST["event_name"]
    details = request.POST["details"]
    severity = request.POST["severity"]

    Event.objects.create(user=request.user,
                         app_name=app_name,
                         name=event_name,
                         date=timezone.now(),
                         details=details,
                         severity=severity)
    return HttpResponse()
