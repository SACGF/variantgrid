from typing import Optional
from dateutil.tz import gettz
from threadlocals.threadlocals import get_current_user, set_request_variable, get_request_variable
from snpdb.models import UserSettings


class UserSettingsManager:

    @staticmethod
    def get_user_settings() -> UserSettings:
        user_settings: Optional[UserSettings]
        if user_settings := get_request_variable("user_settings"):
            return user_settings

        if user := get_current_user():
            user_settings = UserSettings.get_for_user(user)
            set_request_variable("user_settings", user_settings)

        return user_settings

    @staticmethod
    def get_user_timezone():
        if user_settings := UserSettingsManager.get_user_settings():
            return user_settings.tz
        else:
            from django.conf import settings
            gettz(settings.TIME_ZONE)
