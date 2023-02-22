from typing import Optional

from dateutil.tz import gettz
from django.contrib.auth.models import User
from threadlocals.threadlocals import get_current_user, set_request_variable, get_request_variable

from snpdb.models import UserSettings


class UserSettingsManager:

    @staticmethod
    def get_user_settings(user: Optional[User] = None) -> UserSettings:
        """
        :param user: Preferred if this isn't provided (in which case it'll default tothe logged in user)
        :return:
        """
        current_user = get_current_user()
        if user and user != current_user:
            # if we're looking at non-current user, just return now, don't look at cache
            return UserSettings.get_for(user)

        # otherwise, see if we've asked for this before and returned cached value
        user = current_user
        if user_settings := get_request_variable("user_settings"):
            return user_settings

        user_settings = UserSettings.get_for(user)
        set_request_variable("user_settings", user_settings)

        return user_settings

    @staticmethod
    def get_user_timezone():
        if user_settings := UserSettingsManager.get_user_settings():
            return user_settings.tz
        else:
            from django.conf import settings
            gettz(settings.TIME_ZONE)
