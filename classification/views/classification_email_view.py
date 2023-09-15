import collections
import datetime
import unicodedata
from datetime import timezone
from functools import cached_property
from typing import List, Optional

import celery
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.http.request import HttpRequest
from django.http.response import HttpResponse
from django.template.loader import render_to_string
from django.utils.timesince import timesince

from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import Classification, classification_flag_types, \
    DiscordanceReportClassification, DiscordanceReport
from classification.models.discordance_models_utils import DiscordanceReportCategories
from email_manager.models import EmailLog
from flags.models import FlagCollection, Flag
from library.log_utils import report_exc_info, report_message
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab, UserSettings, GenomeBuild

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

    def _get_discordance_report_summaries(self):
        discordant_vcs = FlagCollection.filter_for_open_flags(
            Classification.objects.filter(lab=self.lab),
            flag_types=[classification_flag_types.discordant]
            # if flag types was not set it will be all flag types
        )

        report_ids = DiscordanceReportClassification.objects.filter(
            classification_original__classification__in=discordant_vcs,
            report__resolution=DiscordanceReportResolution.ONGOING).values_list('report', flat=True)
        dr_qs = DiscordanceReport.objects.filter(pk__in=report_ids).order_by('-id')
        return dr_qs

    @cached_property
    def discordance_report_categories(self) -> DiscordanceReportCategories:
        return DiscordanceReportCategories(perspective=LabPickerData.for_lab(self.lab))

    @cached_property
    def flagged_variants(self) -> QuerySet[Flag]:
        return FlagCollection.filter_for_open_flags(
            Classification.objects.filter(lab=self.lab).exclude(withdrawn=True)
        ).order_by('-created')

    @cached_property
    def pending_changes(self) -> QuerySet[Flag]:
        return FlagCollection.filter_for_open_flags(
            Classification.objects.filter(lab=self.lab).exclude(withdrawn=True),
            flag_types=[classification_flag_types.classification_pending_changes]
        ).order_by('-created')

    @property
    def flagged_variants_count(self) -> int:
        return self.flagged_variants.count()

    @property
    def pending_changes_count(self) -> int:
        return self.pending_changes.count()

    @cached_property
    def imported_30_days_count(self):
        since = timezone.now() - datetime.timedelta(days=30)
        vcgs = Classification.objects.filter(lab=self.lab).filter(created__gte=since)
        return vcgs.count()


class EmailSummaryData:

    def __init__(self, perspective: LabPickerData):
        labs = sorted(perspective.selected_labs)
        self.lab_summaries: List[EmailLabSummaryData] = [EmailLabSummaryData(lab=lab, user=perspective.user) for lab in labs]


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


def send_summary_email_to_user(user: User):
    discordance_email = settings.DISCORDANCE_EMAIL
    if discordance_email:
        content = summary_email_content(LabPickerData.for_user(user))

        return EmailLog.send_mail(subject=content.subject,
                                  html=content.html,
                                  text=content.text,
                                  from_email=discordance_email,
                                  recipient_list=[user.email])


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
