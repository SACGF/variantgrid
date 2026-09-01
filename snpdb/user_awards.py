"""
Award definitions registry - see claude/plans/1819_user_flair_and_awards_plan.md

Each app declares its awards in <app>/user_awards.py and imports that module from AppConfig.ready()
(the same way signal receivers are loaded). The nightly/hourly recompute is in
snpdb/user_award_updates.py; this module deliberately imports no models so the UserAward model can
look its definition up.
"""
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from django.conf import settings

from snpdb.models.models_enums import AwardPeriod, UserAwardKind

# {subject: {user_id: count}} - subject None for un-subjected awards
AwardCounts = dict[Optional[str], dict[int, int]]
AwardCounter = Callable[[Optional[datetime]], AwardCounts]


@dataclass(frozen=True)
class AwardDefinition:
    key: str  # "top_tagger", "cold_case"
    kind: UserAwardKind
    title: str  # "Top tagger"
    description: str  # how to earn it - shown on locked badges
    icon: str  # font-awesome class e.g. "fa-tag" - profile, and the subject icon on grids
    counter: AwardCounter
    periods: tuple[AwardPeriod, ...] = ()  # titles only
    tiers: tuple[int, ...] = ()  # badges only: (bronze, silver, gold) thresholds
    hidden: bool = False  # badges only: listed on the profile only once earned

    def __post_init__(self):
        if self.kind == UserAwardKind.TITLE and not self.periods:
            raise ValueError(f"Title award '{self.key}' needs periods")
        if self.kind == UserAwardKind.BADGE and len(self.tiers) != 3:
            raise ValueError(f"Badge award '{self.key}' needs (bronze, silver, gold) tiers")

    def compute(self, since: Optional[datetime]) -> AwardCounts:
        """ since=None means all time. Badges ignore since """
        return self.counter(since)

    @property
    def enabled(self) -> bool:
        return self.key not in settings.USER_AWARDS_DISABLED_KEYS

    def award_text(self, subject: Optional[str] = None, period: Optional[AwardPeriod] = None) -> str:
        text = self.title
        if subject:
            text += f" of '{subject}'"
        if period:
            text += f" ({AwardPeriod(period).label})"
        return text


_AWARD_DEFINITIONS: dict[str, AwardDefinition] = {}


def register_award(definition: AwardDefinition) -> AwardDefinition:
    _AWARD_DEFINITIONS[definition.key] = definition
    return definition


def get_award_definitions() -> list[AwardDefinition]:
    return list(_AWARD_DEFINITIONS.values())


def get_award_definition(key: str) -> Optional[AwardDefinition]:
    return _AWARD_DEFINITIONS.get(key)


def get_badge_definitions() -> list[AwardDefinition]:
    return [d for d in get_award_definitions() if d.kind == UserAwardKind.BADGE]


def get_title_definitions() -> list[AwardDefinition]:
    return [d for d in get_award_definitions() if d.kind == UserAwardKind.TITLE]


def counts_by_user(rows, subject_key: Optional[str] = None) -> AwardCounts:
    """ Turn a .values(user_key, [subject_key]).annotate(count=...) queryset into AwardCounts.
        Rows need 'user_id' and 'count' keys. With subject_key, each row also lands in the
        overall (None) subject so 'top overall' comes for free alongside 'top per subject' """
    counts: AwardCounts = {}
    for row in rows:
        user_id = row["user_id"]
        count = row["count"]
        if subject_key:
            counts.setdefault(row[subject_key], {})[user_id] = count
            overall = counts.setdefault(None, {})
            overall[user_id] = overall.get(user_id, 0) + count
        else:
            counts.setdefault(None, {})[user_id] = count
    return counts
