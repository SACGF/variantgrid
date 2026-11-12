"""
API tokens: ApiToken is the knox token model (settings.KNOX_TOKEN_MODEL) - several per user, each with a name and an
expiry, with only a digest of the key stored - and ApiTokenUse is its use history, one row per hour a token is used.
Create one with `ApiToken.create_for_user` (the key is returned once, never stored) and end it with `ApiToken.revoke`;
both write an audit Event. snpdb/api_token_authentication.py:ApiTokenAuthentication authenticates them and calls
`ApiToken.record_use`.
"""
from datetime import timedelta
from typing import Optional

from django.contrib.auth.models import User
from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE, SET_NULL
from django.utils import timezone
from knox.models import AbstractAuthToken

from eventlog.models import create_event


class ApiToken(AbstractAuthToken):
    """ Expired and revoked tokens are kept, so their use history stays for an admin to review """
    USE_RECORD_INTERVAL = timedelta(hours=1)

    name = models.TextField()
    last_used = models.DateTimeField(null=True, blank=True)
    revoked = models.DateTimeField(null=True, blank=True)
    revoked_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL, related_name='+')

    def __str__(self):
        return f"{self.user}: {self.name}"

    @staticmethod
    def create_for_user(user: User, name: str, expiry: Optional[timedelta],
                        created_by: Optional[User]) -> tuple['ApiToken', str]:
        """ Returns the token and its key - the key is not stored, so this is the only time it is available """
        token, key = ApiToken.objects.create(user=user, expiry=expiry, name=name)
        expires = f"expires {token.expiry:%Y-%m-%d}" if token.expiry else "never expires"
        create_event(created_by, "api_token_created", details=f"Token '{name}' for {user.username}, {expires}",
                     app_name="snpdb")
        return token, key

    @property
    def is_expired(self) -> bool:
        return self.expiry is not None and self.expiry <= timezone.now()

    @property
    def is_valid(self) -> bool:
        return self.revoked is None and not self.is_expired

    @property
    def status(self) -> str:
        if self.revoked:
            return "Revoked"
        if self.is_expired:
            return "Expired"
        return "Active"

    def can_manage(self, user: User) -> bool:
        return user == self.user or user.is_superuser

    def revoke(self, user: User):
        if self.revoked:
            return
        self.revoked = timezone.now()
        self.revoked_by = user
        self.save(update_fields=["revoked", "revoked_by"])
        create_event(user, "api_token_revoked", details=f"Token '{self.name}' for {self.user.username}",
                     app_name="snpdb")

    def record_use(self, ip_address: Optional[str], path: str):
        """ Records at most one use per USE_RECORD_INTERVAL: the conditional update means requests arriving
            together write a single ApiTokenUse """
        now = timezone.now()
        last_recorded_before = now - self.USE_RECORD_INTERVAL
        if self.last_used and self.last_used > last_recorded_before:
            return
        not_recently_used = Q(last_used__isnull=True) | Q(last_used__lte=last_recorded_before)
        if ApiToken.objects.filter(not_recently_used, pk=self.pk).update(last_used=now):
            self.last_used = now
            ApiTokenUse.objects.create(token=self, used=now, ip_address=ip_address, path=path)


class ApiTokenUse(models.Model):
    """ The first request a token made in each ApiToken.USE_RECORD_INTERVAL """
    token = models.ForeignKey(ApiToken, on_delete=CASCADE, related_name='uses', db_index=False)
    used = models.DateTimeField()
    ip_address = models.GenericIPAddressField(null=True)
    path = models.TextField()

    class Meta:
        indexes = [models.Index(fields=["token", "used"])]

    def __str__(self):
        return f"{self.token}: {self.used:%Y-%m-%d %H:%M} {self.path}"
