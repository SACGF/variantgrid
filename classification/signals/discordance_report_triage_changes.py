import logging

from django.db.models.signals import post_save
from django.dispatch import receiver

from classification.models import DiscordanceReportTriageStatus, DiscordanceReport, discordance_change_signal, \
    ensure_discordance_report_triages_for, DiscordanceReportTriage


"""
Creates, closes or re-opens triages base don activity of the discordance
"""


@receiver(discordance_change_signal, sender=DiscordanceReport)
def tidy_triage(sender, discordance_report: DiscordanceReport, **kwargs):
    ensure_discordance_report_triages_for(discordance_report)
