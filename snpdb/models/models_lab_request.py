"""
Labs for self-registered users, and requests to turn them into real ones.

settings.USER_CREATE_ORG_LABS gives each new user a lab in an organisation (create_user_org_labs, run on user
creation). A lab pattern filled from the user ("user_group_%(username)s") makes it a private lab in a holding
organisation: users join their own lab's group but not the organisation's, so nothing they upload or share with the
organisation reaches other self-registered users.

A LabRequest asks the admins to make the requester's private lab their institution's lab. The requester proves they
own an email on the institution's website domain (email_matches_website, then the link from
send_confirmation_email), which notifies the admins. An admin checks the institution by hand and approves
(snpdb/lab_request_approval.py:approve_lab_request, which renames the private lab in place) or rejects it.
Views are in snpdb/views/views_lab_request.py.
"""
from datetime import timedelta
from typing import Optional
from urllib.parse import urlparse

from django.conf import settings
from django.contrib.auth.models import User
from django.core import signing
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL
from django.urls import reverse
from django.utils import timezone
from django.utils.text import slugify
from model_utils.models import TimeStampedModel

from email_manager.models.email_model import EmailLog
from eventlog.notifications import AdminNotificationBuilder
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from snpdb.enums import LabRequestStatus
from snpdb.models.models import Lab, Organization
from user_messages.models import Message


def is_per_user_lab_pattern(lab_pattern: str) -> bool:
    return "%(" in lab_pattern


def _user_lab_name(lab_pattern: str, user: User) -> str:
    """ Fields read through getattr, as request.user is a lazy wrapper whose __dict__ isn't the user's """
    user_values = {field.attname: getattr(user, field.attname) for field in User._meta.concrete_fields}
    return lab_pattern % user_values


def user_create_org_lab(org_group_name: str, lab_pattern: str, user: User) -> Lab:
    """ The settings.USER_CREATE_ORG_LABS lab for user (and its organisation), created if missing. Doesn't add
        user to it """
    organization, _ = Organization.objects.get_or_create(group_name=org_group_name)
    lab_name = _user_lab_name(lab_pattern, user)
    lab, _ = Lab.objects.get_or_create(group_name=f"{organization.group_name}/{lab_name.lower()}",
                                       defaults={"name": lab_name, "organization": organization})
    return lab


def create_user_org_labs(user: User):
    """ Adds a new user to each settings.USER_CREATE_ORG_LABS lab. A private lab's user doesn't join the holding
        organisation's group """
    for org_group_name, lab_pattern in settings.USER_CREATE_ORG_LABS.items():
        lab = user_create_org_lab(org_group_name, lab_pattern, user)
        lab.group.user_set.add(user)
        if is_per_user_lab_pattern(lab_pattern):
            lab.organization.send_user_create_message(user)
        else:
            lab.organization.add_member(user)


def get_user_private_lab(user: User) -> Optional[Lab]:
    """ The user's private lab while it is still in its holding organisation - approving a LabRequest moves it out """
    for org_group_name, lab_pattern in settings.USER_CREATE_ORG_LABS.items():
        if is_per_user_lab_pattern(lab_pattern):
            private_group_name = f"{org_group_name}/{_user_lab_name(lab_pattern, user).lower()}"
            if lab := Lab.objects.filter(group_name=private_group_name, organization__group_name=org_group_name).first():
                return lab
    return None


def group_name_slug(name: str) -> str:
    return slugify(name).replace("-", "_")


def lab_group_name(organization_group_name: str, lab_name: str) -> str:
    return f"{organization_group_name}/{group_name_slug(lab_name)}"


def email_matches_website(email: str, website: str) -> bool:
    """ The email's domain is the website's host (ignoring www.), or one is a subdomain of the other:
        someone@rsu.lv matches www.rsu.lv and genetics.rsu.lv, someone@mail.rsu.lv matches rsu.lv """
    email_domain = email.rsplit("@", 1)[-1].lower()
    host = (urlparse(website).hostname or "").lower().removeprefix("www.")
    if not email_domain or not host:
        return False
    return email_domain == host or email_domain.endswith(f".{host}") or host.endswith(f".{email_domain}")


class LabRequest(TimeStampedModel):
    CONFIRMATION_SALT = "snpdb.lab_request.confirm"
    CONFIRMATION_MAX_AGE = timedelta(days=7)

    user = models.ForeignKey(User, on_delete=CASCADE)
    organization_name = models.TextField()
    lab_name = models.TextField()
    website = models.URLField()
    email = models.EmailField()
    email_confirmed = models.DateTimeField(null=True, blank=True)
    status = models.CharField(max_length=1, choices=LabRequestStatus.choices, default=LabRequestStatus.PENDING)
    reviewed_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL, related_name="+")
    reviewed = models.DateTimeField(null=True, blank=True)
    lab = models.ForeignKey(Lab, null=True, blank=True, on_delete=SET_NULL)  # Set on approval

    class Meta:
        ordering = ["-created"]

    def __str__(self):
        return f"{self.user}: {self.lab_name} ({self.organization_name})"

    def get_absolute_url(self):
        return reverse("view_lab_request", kwargs={"pk": self.pk})

    @property
    def can_review(self) -> bool:
        return self.status == LabRequestStatus.PENDING and self.email_confirmed is not None

    @property
    def confirmation_token(self) -> str:
        return signing.dumps({"pk": self.pk, "email": self.email}, salt=self.CONFIRMATION_SALT)

    @staticmethod
    def from_confirmation_token(token: str) -> 'LabRequest':
        """ Raises signing.BadSignature (or its subclass SignatureExpired) for a bad or expired token """
        data = signing.loads(token, salt=LabRequest.CONFIRMATION_SALT,
                             max_age=LabRequest.CONFIRMATION_MAX_AGE.total_seconds())
        try:
            return LabRequest.objects.get(pk=data["pk"], email=data["email"])
        except LabRequest.DoesNotExist as e:
            raise signing.BadSignature("Lab request no longer exists") from e

    def send_confirmation_email(self):
        confirm_url = get_url_from_view_path(reverse("lab_request_confirm", kwargs={"token": self.confirmation_token}))
        # No text from the request: anyone can enter any address, so this email carries nothing they wrote
        text = (f"A {settings.SITE_NAME} user asked to set up a lab using this email address.\n\n"
                f"To confirm this is your address, open this link (it expires in {self.CONFIRMATION_MAX_AGE.days} days):\n"
                f"{confirm_url}\n\nIf you didn't ask for this, ignore this email.")
        EmailLog.send_mail(subject=f"Confirm your {settings.SITE_NAME} lab request",
                           html="", text=text, from_email=settings.DEFAULT_FROM_EMAIL, recipient_list=[self.email])

    def confirm_email(self):
        """ Records the requester owns the email, and asks the admins to review the request """
        if self.email_confirmed:
            return
        self.email_confirmed = timezone.now()
        self.save(update_fields=["email_confirmed", "modified"])

        # The request's text is the user's, so it is shown only on the review page, where templates escape it
        nb = AdminNotificationBuilder(message="Lab request", is_communication=True)
        nb.add_markdown(f"{self.user.username} confirmed their email for a lab request. "
                        f"[Review the request]({get_url_from_view_path(self.get_absolute_url())})")
        nb.send()

    def set_reviewed(self, status: LabRequestStatus, reviewed_by: User):
        self.status = status
        self.reviewed_by = reviewed_by
        self.reviewed = timezone.now()

    def notify_requester(self, subject: str, body: str):
        Message.objects.create(subject=subject, body=body, sender=admin_bot(), recipient=self.user)
        EmailLog.send_mail(subject=subject, html="", text=body,
                           from_email=settings.DEFAULT_FROM_EMAIL, recipient_list=[self.email])

    def reject(self, reviewed_by: User, notify: bool):
        self.set_reviewed(LabRequestStatus.REJECTED, reviewed_by)
        self.save()
        if notify:
            self.notify_requester(f"{settings.SITE_NAME} lab request",
                                  f"Your request for the lab '{self.lab_name}' at '{self.organization_name}' "
                                  f"was not approved.")
