"""
Award definitions registry

Each app declares its badges in <app>/user_awards.py and imports that module from AppConfig.ready()
(the same way signal receivers are loaded). The nightly recompute is in snpdb/user_award_updates.py;
this module deliberately imports no models so the UserAward model can look its definition up.
"""
from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings

# {user_id: count}
AwardCounts = dict[int, int]
AwardCounter = Callable[[], AwardCounts]


@dataclass(frozen=True)
class AwardDefinition:
    key: str  # "tagger", "cold_case"
    title: str  # "Tagger"
    description: str  # how to earn it - shown on locked badges
    icon: str  # font-awesome class e.g. "fa-tag"
    counter: AwardCounter
    tiers: tuple[int, ...]  # (bronze, silver, gold) thresholds

    def __post_init__(self):
        if len(self.tiers) != 3:
            raise ValueError(f"Badge award '{self.key}' needs (bronze, silver, gold) tiers")

    def compute(self) -> AwardCounts:
        return self.counter()

    @property
    def enabled(self) -> bool:
        return self.key not in settings.USER_AWARDS_DISABLED_KEYS


_AWARD_DEFINITIONS: dict[str, AwardDefinition] = {}


def register_award(definition: AwardDefinition) -> AwardDefinition:
    _AWARD_DEFINITIONS[definition.key] = definition
    return definition


def get_award_definitions() -> list[AwardDefinition]:
    return list(_AWARD_DEFINITIONS.values())


def get_award_definition(key: str) -> Optional[AwardDefinition]:
    return _AWARD_DEFINITIONS.get(key)


def counts_by_user(rows) -> AwardCounts:
    """ Turn a .values("user_id").annotate(count=...) queryset into AwardCounts """
    return {row["user_id"]: row["count"] for row in rows}
