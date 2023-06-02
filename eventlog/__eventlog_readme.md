# Event Log

The Event Log is invoked to record some actions (not used as much as it once was).

This also implements ViewEvents, which are automatically recorded when a URL is accessed (see middleware.py).
To record ViewEvents, settings.LOG_ACTIVITY_APPS does have to be set.