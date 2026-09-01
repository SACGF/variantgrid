"""
snpdb's own award definitions (#1819) - registered from SnpdbConfig.ready()
"""
from datetime import datetime
from typing import Optional

from django.contrib.auth.models import User
from django.utils import timezone

from snpdb.models.models_enums import UserAwardKind
from snpdb.user_awards import AwardCounts, AwardDefinition, register_award


def _years_since_joined(_since: Optional[datetime]) -> AwardCounts:
    now = timezone.now()
    counts = {}
    for user_id, date_joined in User.objects.filter(is_active=True).values_list("id", "date_joined"):
        if years := (now - date_joined).days // 365:
            counts[user_id] = years
    return {None: counts}


register_award(AwardDefinition(
    key="founder",
    kind=UserAwardKind.BADGE,
    title="Founder",
    description="Years as a member",
    icon="fa-landmark",
    counter=_years_since_joined,
    tiers=(1, 3, 5),
))
