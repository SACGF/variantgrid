"""
Wiki award definitions (#1819) - gene, gene list, variant, node and sequencing run wikis all share
the snpdb Wiki base. Registered from GenesConfig.ready()
"""
from django.db.models import Count, F

from snpdb.models import Wiki
from snpdb.user_awards import AwardCounts, AwardDefinition, counts_by_user, register_award


def _wikis_edited() -> AwardCounts:
    """ Wikis have no edit history, so this counts wikis whose last edit was by the user - an
        approximation that undercounts busy pages """
    qs = Wiki.objects.filter(last_edited_by__isnull=False)
    return counts_by_user(qs.values(user_id=F("last_edited_by")).annotate(count=Count("id")))


register_award(AwardDefinition(
    key="wiki_scribe",
    title="Wiki scribe",
    description="Wiki pages written",
    icon="fa-feather",
    counter=_wikis_edited,
    tiers=(5, 25, 100),
))
