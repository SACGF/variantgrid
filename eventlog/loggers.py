from django.utils import timezone
import logging
import re

from library.enums.log_level import LogLevel
from library.utils import import_class


missing_mapping_re = re.compile(r'Not Found:.*\.js.*\.map')

class EventLogHandler(logging.Handler):  # Inherit from logging.Handler

    def __init__(self):
        logging.Handler.__init__(self)

    def emit(self, record):
        try:
            # Have to do it this way as logging is loaded in settings.py, and if you import a model then you get:
            # Model class eventlog.models.Event doesn't declare an explicit app_label and either isn't in an application in INSTALLED_APPS
            # or else was imported before its application was loaded. This will no longer be supported in Django 1.9.
            Event = import_class('eventlog.models.Event')

            user = record.request.user
            app_name = record.module.split('.')[0]
            if record.levelname == 'WARNING':
                severity = LogLevel.WARNING
            elif record.levelname == 'INFO':
                severity = LogLevel.INFO
            elif record.levelname == 'ERROR':
                severity = LogLevel.ERROR

            name = 'django_exception'
            details_list = []
            for value in [record.getMessage(), record.exc_text]:
                if value:
                    details_list.append(value)

            details = '\n'.join(map(str, details_list))
            # don't fill up the logs with missing map files
            if missing_mapping_re.match(details):
                return

            Event.objects.create(user=user,
                                 app_name=app_name,
                                 name=name,
                                 date=timezone.now(),
                                 details=details,
                                 severity=severity)
            #rollbar.init(ROLLBAR)
            #rollbar.report_message(message = name,
            #                       level = severity.lower(),
            #                       request = record.request,
            #                       extra_data = {'name': name, 'app_name': app_name})
        except:
            pass

        return
