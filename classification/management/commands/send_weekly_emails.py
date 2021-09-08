from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from classification.views.classification_email_view import send_summary_email_to_user
from library.log_utils import report_exc_info
from snpdb.models import UserSettings


class Command(BaseCommand):

    def handle(self, *args, **options):
        for user in User.objects.filter(is_active=True):
            try:
                us = UserSettings.get_for_user(user)
                if us.get_for_user(user).email_weekly_updates:
                    send_summary_email_to_user(user=user)
            except:
                report_exc_info({"user": user.username})
