"""
LogLevel: the severity strings (DEBUG .. ERROR) stored on eventlog Event rows and used by
notification builders; `vg status` filters ERROR.
"""


class LogLevel:
    DEBUG = 'D'
    INFO = 'I'
    WARNING = 'W'
    ERROR = 'E'

    CHOICES = [(DEBUG, 'DEBUG'),
               (INFO, 'INFO'),
               (WARNING, 'WARNING'),
               (ERROR, 'ERROR')]
