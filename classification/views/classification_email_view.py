import collections
import unicodedata
from functools import cached_property
from typing import Optional
import celery
from django.conf import settings
from django.contrib.auth.models import User
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.template.loader import render_to_string
from django.utils.timesince import timesince
from classification.models import Classification
from classification.services.overlaps_services import OverlapsSummary
from email_manager.models import EmailLog
from library.log_utils import report_exc_info, report_message
from snpdb.lab_picker import LabPickerData
from snpdb.models import GenomeBuild, Lab, UserSettings

EmailOutput = collections.namedtuple('EmailOutput', 'subject html text')


class EmailLabSummaryData:

    def __init__(self, lab: Lab, user: User):
        self.lab = lab
        self.user = user

    @cached_property
    def last_imported_new_ago(self) -> Optional[str]:
        latest = Classification.objects.order_by('-created').filter(lab=self.lab).values_list('created',
                                                                                              flat=True).first()
        if latest:
            timesince_str = timesince(latest)
            return unicodedata.normalize("NFKD", timesince_str)
        else:
            return None

    @cached_property
    def genome_build(self) -> GenomeBuild:
        # TODO if user setting isn't set, grab it from lab instead
        return UserSettings.get_genome_build_or_default(self.user)

    @cached_property
    def overlaps_summary(self):
        return OverlapsSummary(perspective=LabPickerData.for_lab(self.lab))


class EmailSummaryData:

    def __init__(self, perspective: LabPickerData):
        labs = sorted(perspective.selected_labs)
        self.lab_summaries: list[EmailLabSummaryData] = [EmailLabSummaryData(lab=lab, user=perspective.user) for lab in labs]


@celery.shared_task
def send_summary_emails():
    report_message("Attempting to send weekly summary emails", level="info")
    for user in User.objects.filter(is_active=True):
        try:
            us = UserSettings.get_for_user(user)
            if us.email_weekly_updates:
                send_summary_email_to_user(user=user)
        except Exception:
            report_exc_info({"user": user.username})


def send_summary_email_to_user(user: User) -> bool:
    discordance_email = settings.DISCORDANCE_EMAIL
    if discordance_email:
        if not Lab.has_active_lab(user, admin_check=True):
            # User has no accessible lab, so there's nothing to summarise - skip rather than raise
            report_message("Skipping weekly summary email - user has no accessible lab", level="debug",
                           extra_data={"user": user.username})
            return False

        content = summary_email_content(LabPickerData.for_user(user))

        return EmailLog.send_mail(subject=content.subject,
                                  html=content.html,
                                  text=content.text,
                                  from_email=discordance_email,
                                  recipient_list=[user.email])
    else:
        return False


def summary_email_preview_html(request: HttpRequest, lab_id: Optional[str] = None) -> HttpResponse:
    return HttpResponse(
        summary_email_content(LabPickerData.for_user(request.user, selection=lab_id)).html
    )


def summary_email_preview_text(request: HttpRequest, lab_id: Optional[str] = None) -> HttpResponse:
    return HttpResponse(
        summary_email_content(LabPickerData.for_user(request.user, selection=lab_id)).text,
        content_type="text/plain"
    )


def summary_email_content(perspective: LabPickerData) -> EmailOutput:
    data = EmailSummaryData(perspective=perspective)
    subject = 'Weekly Classification Summary'

    context = {
        "data": data
    }

    html = render_to_string('classification/emails/classification_summary_email.html', context).strip()
    text = render_to_string('classification/emails/classification_summary_email.txt', context).strip()
    return EmailOutput(subject=subject, html=html, text=text)
