"""
eventlog's enums as one namespace, outside eventlog/models.py so importing one doesn't load the models.
Callers write `from eventlog.enums import IntegrationActivityStatus`.
"""
from eventlog.enums.eventlog_enums import *
