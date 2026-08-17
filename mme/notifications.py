""" Notify the depositor - the curator who owns a matched classification - as well as the
querier, for both the inbound and outbound halves of MatchMaker Exchange.

Peers re-query on a schedule, so the same (our classification, their patient) pair recurs;
we notify once per pair. Deliberately not routed through LabNotificationBuilder
(`snpdb/utils.py`), which hardcodes `settings.DISCORDANCE_EMAIL` as its from-address - only
set on Shariant, so elsewhere it silently sends nothing. Its recipient fallback is worth
copying, and `lab_notification_emails` does.
"""
import logging

from django.conf import settings
from django.contrib.auth.models import Group
from django.urls import reverse
from django.utils import timezone

from email_manager.models import EmailLog
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from library.log_utils import AdminNotificationBuilder
from mme.models import MMEInboundMatch
from user_messages.models import Message


def lab_notification_emails(lab) -> list[str]:
    """ `lab.email` is the lab's notification inbox; failing that every active member, so a
        match still lands when it is unset or the authoring curator has gone. """
    if lab is None:
        return []
    if lab_email := (lab.email or "").strip():
        return [lab_email]
    if not lab.group_name:
        return []
    return list(Group.objects.get(name=lab.group_name)
                .user_set.filter(is_active=True)
                .exclude(email="")
                .exclude(email__isnull=True)
                .values_list("email", flat=True))


def notify_match(classification, subject: str, body: str) -> bool:
    """ Notify the authoring curator AND the lab - a case can sit unmatched for years, so
        the lab is the durable recipient. True if any channel accepted it; never raises. """
    sent = False

    user = classification.user
    if user is not None and user.is_active:     # inactive = left; the lab still hears
        try:
            Message.objects.create(subject=subject, body=body,
                                   sender=admin_bot(), recipient=user)
            sent = True
        except Exception:
            logging.exception("MME: in-app notification failed for classification %s",
                              classification.pk)

    if emails := lab_notification_emails(getattr(classification, "lab", None)):
        try:
            EmailLog.send_mail(subject=subject, html="", text=body,
                               from_email=settings.MME_FROM_EMAIL, recipient_list=emails)
            sent = True
        except Exception:
            logging.exception("MME: lab email failed for classification %s", classification.pk)

    if not sent:
        # A match nobody can be told about is an incident, not a log line.
        nb = AdminNotificationBuilder("MME match notification failed")
        nb.add_markdown(f"No notification channel accepted the MatchMaker Exchange match "
                        f"notification for classification {classification.pk}. Check "
                        f"`MME_FROM_EMAIL`, `SEND_EMAILS`, and the lab's notification email.")
        nb.send()
    return sent


def record_inbound_matches(inbound_query, matches, query_patient: dict) -> list[MMEInboundMatch]:
    """ Persist which of our classifications we handed to the peer (audit + follow-up). """
    remote_patient_id = str(query_patient.get("id") or "")
    return MMEInboundMatch.objects.bulk_create([
        MMEInboundMatch(
            inbound_query=inbound_query,
            classification=match.classification,
            score=match.score,
            remote_patient_id=remote_patient_id,
            query_patient_json=query_patient,
        ) for match in matches
    ])


def _already_notified(match: MMEInboundMatch) -> bool:
    """ True when this curator has already been told about this remote case. """
    return (MMEInboundMatch.objects
            .filter(classification_id=match.classification_id,
                    remote_patient_id=match.remote_patient_id,
                    notified__isnull=False)
            .exclude(pk=match.pk)
            .exists())


def notify_depositors(inbound_query_id: int) -> None:
    """ Tell the depositors of each freshly matched classification. Best-effort: a
        messaging failure must not fail the peer's already-answered /match request. """
    for match in (MMEInboundMatch.objects
                  .filter(inbound_query_id=inbound_query_id, notified__isnull=True)
                  .select_related("classification__lab", "inbound_query")):
        if match.remote_patient_id and _already_notified(match):
            continue
        peer = match.inbound_query.peer_node_id or "an MME peer"
        contact = (match.query_patient_json.get("contact") or {})
        contact_line = ""
        if contact.get("name"):
            href = contact.get("href") or ""
            contact_line = f"\n\nTheir contact: {contact['name']}" + (f" — <{href}>" if href else "")

        url = get_url_from_view_path(
            reverse("mme_view_inbound_match", kwargs={"inbound_match_id": match.pk}))
        sent = notify_match(
            classification=match.classification,
            subject=f"MatchMaker Exchange: {peer} matched your classification",
            body=(f"`{peer}` queried MatchMaker Exchange and one of your classifications "
                  f"(#{match.classification_id}) was returned as a possible match "
                  f"(score {match.score:.2f})."
                  f"{contact_line}\n\n[View the match]({url})."),
        )
        if sent:
            match.notified = timezone.now()
            match.save(update_fields=["notified"])
