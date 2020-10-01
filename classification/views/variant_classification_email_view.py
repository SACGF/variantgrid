import collections
import datetime
from typing import List

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.template.loader import render_to_string
from django.utils.timesince import timesince
from lazy import lazy

from email_manager.models import EmailLog
from flags.models import FlagCollection
from library.log_utils import report_exc_info, report_message
from snpdb.models import Lab, UserSettings
from classification.enums import ShareLevel
from classification.enums.discordance_enums import DiscordanceReportResolution
from classification.models import VariantClassification, variant_classification_flag_types, \
    DiscordanceReportClassification, DiscordanceReport
import celery
from django.http.request import HttpRequest
from django.http.response import HttpResponse


EmailOutput = collections.namedtuple('EmailOutput', 'subject html text')


class EmailLabSummaryData:

    def __init__(self, lab: Lab, user: User):
        self.lab = lab
        self.user = user

    @lazy
    def last_imported_new_ago(self):
        latest = VariantClassification.objects.order_by('-created').filter(lab=self.lab).values_list('created',
                                                                                                flat=True).first()
        if latest:
            return timesince(latest)
        else:
            return None

    @lazy
    def discordance_reports(self) -> List[DiscordanceReport]:
        discordant_vcs = FlagCollection.filter_for_open_flags(
            VariantClassification.objects.filter(lab=self.lab),
            flag_types=[variant_classification_flag_types.discordant]
            # if flag types was not set it will be all flag types
        )

        report_ids = DiscordanceReportClassification.objects.filter(
            classification_original__variant_classification__in=discordant_vcs,
            report__resolution=DiscordanceReportResolution.ONGOING).values_list('report', flat=True)
        return DiscordanceReport.objects.filter(pk__in=report_ids).order_by('id')

    @lazy
    def flagged_variants(self):
        vcqs = FlagCollection.filter_for_open_flags(
            VariantClassification.objects.filter(lab=self.lab).exclude(withdrawn=True)
        ).order_by('-created')

        return vcqs

    @lazy
    def flagged_variants_count(self):
        return self.flagged_variants.count()

    @lazy
    def imported_30_days_count(self):
        since = datetime.datetime.now() - datetime.timedelta(days=30)
        vcgs = VariantClassification.objects.filter(lab=self.lab).filter(created__gte=since)
        return vcgs.count()


class EmailSummaryData:

    def __init__(self, user: User):
        self.user = user
        # labs
        labs_qs: QuerySet
        if user.is_superuser:
            labs_qs = Lab.objects
        else:
            labs_qs = Lab.valid_labs_qs(user)

        labs_qs = labs_qs.order_by('group_name')
        self.lab_summaries: List[EmailLabSummaryData] = [EmailLabSummaryData(lab=lab, user=user) for lab in labs_qs]

    @property
    def has_issues(self) -> bool:
        for lab_summary in self.lab_summaries:
            if lab_summary.flagged_variants_count > 0:
                return True
        return False


@celery.task
def send_summary_emails():
    report_message("Attempting to send weekly summary emails", level="info")
    for user in User.objects.filter(is_active=True):
        try:
            us = UserSettings.get_for_user(user)
            if us.email_weekly_updates:
                send_summary_email_to_user(user=user)
        except:
            report_exc_info({"user": user.username})


def send_summary_email_to_user(user: User):
    discordance_email = settings.DISCORDANCE_EMAIL
    if discordance_email:
        content = summary_email_content(user)

        return EmailLog.send_mail(subject=content.subject,
                                  html=content.html,
                                  text=content.text,
                                  from_email=discordance_email,
                                  recipient_list=[user.email])


def summary_email_preview_html(request: HttpRequest) -> HttpResponse:
    return HttpResponse(summary_email_content(request.user).html)


def summary_email_preview_text(request: HttpRequest) -> HttpResponse:
    return HttpResponse(summary_email_content(request.user).text, content_type="text/plain")


def summary_email_content(user: User) -> EmailOutput:
    data = EmailSummaryData(user=user)
    subject = 'Weekly Classification Summary'

    context = {
        "data": data
    }

    html = render_to_string('classification/emails/variant_classification_summary_email.html', context).strip()
    text = render_to_string('classification/emails/variant_classification_summary_email.txt', context).strip()
    return EmailOutput(subject=subject, html=html, text=text)
