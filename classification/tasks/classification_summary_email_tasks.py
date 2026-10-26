"""
Weekly classification summary email, sent by the Monday `discordance-emails-weekly` celery beat entry when
DISCORDANCE_EMAIL is set. The email itself is built in `classification.views.classification_email_view`, which also
serves its previews.

Entry point: `send_summary_emails` (celery task).
"""
import celery
from django.contrib.auth.models import User

from classification.views.classification_email_view import send_summary_email_to_user
from library.log_utils import report_exc_info, report_message
from snpdb.models import UserSettings


@celery.shared_task(queue='db_workers')
def send_summary_emails():
    report_message("Attempting to send weekly summary emails", level="info")
    for user in User.objects.filter(is_active=True):
        try:
            us = UserSettings.get_for_user(user)
            if us.email_weekly_updates:
                send_summary_email_to_user(user=user)
        except Exception:
            report_exc_info({"user": user.username})
