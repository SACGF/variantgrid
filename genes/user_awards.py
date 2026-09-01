"""
Wiki award definitions (#1819) - gene, gene list, variant, node and sequencing run wikis all share
the snpdb Wiki base. Registered from GenesConfig.ready()
"""
from datetime import datetime
from typing import Optional

from django.db.models import Count, F

from snpdb.models import Wiki
from snpdb.models.models_enums import AwardPeriod, UserAwardKind
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award


def _wikis_edited(since: Optional[datetime]) -> AwardCounts:
    """ Wikis have no edit history, so this counts wikis whose last edit was by the user in the
        period - an approximation that undercounts busy pages """
    qs = Wiki.objects.filter(last_edited_by__isnull=False)
    if since:
        qs = qs.filter(modified__gte=since)
    return counts_by_user(qs.values(user_id=F("last_edited_by")).annotate(count=Count("id")))


register_award(AwardDefinition(
    key="top_wiki",
    kind=UserAwardKind.TITLE,
    title="Top wiki editor",
    description="Most wiki pages edited",
    icon="fa-book",
    counter=_wikis_edited,
    periods=(AwardPeriod.ALL_TIME, AwardPeriod.MONTH, AwardPeriod.DAY),
))

register_award(AwardDefinition(
    key="wiki_scribe",
    kind=UserAwardKind.BADGE,
    title="Wiki scribe",
    description="Wiki pages written",
    icon="fa-feather",
    counter=_wikis_edited,
    tiers=(5, 25, 100),
))
